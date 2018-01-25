; This is a 3D example of ctl file for DOS-calculation programs. For more information, please refer to our website:
; https://github.com/boyuanliuoptics/DOS-calculation/edit/master/DOS_GGR.m

;-----------------------------------User changing part begins-----------------------------------------

;-------------------
; structure setting
;-------------------

; example of a double gyroid photonic crystal

;set the primitive lattice
(set! geometry-lattice (make lattice
(basis-size (/ (sqrt 3) 2) (/ (sqrt 3) 2) (/ (sqrt 3) 2))
(basis1 -1 1 1)
(basis2 1 -1 1)
(basis3 1 1 -1)))

;define a couple of parameters about the structure (which we can set from the command-line)
(define-param eps 16)
(define-param epsbg 1); the dielectric of air
(define-param isoval -1.1);isovalue to choose packing fraction
(define pi (* 4 (atan 1)))
(define FB
(lambda (h k l x y z)
(* (cos (* 2 pi (+ h k l) 0.25)) (+ 
	(* (sin (* 2 pi (+ (* h x) (/ l 4)))) (sin (* 2 pi (+ (* k y) (/ h 4)))) (sin (* 2 pi (+ (* l z) (/ k 4)))))
	(* (sin (* 2 pi (+ (* h y) (/ l 4)))) (sin (* 2 pi (+ (* k z) (/ h 4)))) (sin (* 2 pi (+ (* l x) (/ k 4)))))
	(* (sin (* 2 pi (+ (* h z) (/ l 4)))) (sin (* 2 pi (+ (* k x) (/ h 4)))) (sin (* 2 pi (+ (* l y) (/ k 4)))))
))))
(define-param b 1)
(define-param c 1)
(define (eps-func p-lattice)
(let* ((p (lattice->cartesian p-lattice))
(x (vector3-x p))
(y (vector3-y p))
(z (vector3-z p)))
(if 
(< (FB 1 1 0 x (/ y b) (/ z c)) isoval)
	(make dielectric-anisotropic (epsilon-diag eps eps eps) (epsilon-offdiag 0 0 0))
	(if (< (FB 1 1 0 (- x) (- (/ y b)) (- (/ z c))) isoval)
		(make dielectric-anisotropic (epsilon-diag eps eps eps) (epsilon-offdiag 0 0 0))
		(make dielectric (epsilon epsbg)))
)
))

(set! default-material (make material-function (material-func eps-func)))

;-----------------------------------------------------
; computing setting -- for drawing band stucture line
;-----------------------------------------------------

;switch for computing the band in the line 
;if you want only want to calculate DOS, change 'true' to any other word. It will be detected in the shell file.
;switch-> true

;high symmetry points in the band
(define G (vector3 0 0 0))
(define H (vector3 0.5 -0.5 0.5))
(define P (vector3 0.75 -0.25 -0.25))
(define N (vector3 0.5 0. -0.5))

;sequence of band structure in line
(define bandseq (list H G N P G))

;k-int is the inter number of k points in line between two high symmetry points
(define-param k-int 30)

;--------------------------------------------------------
; computing setting -- for density of states computation
;--------------------------------------------------------

;when there is time-reversial (T) symmetry we can use half of the Brillouin zone (BZ), that is x ranges in [0,bx/2] with num-k/2-1 points in this dimension
;the 'true' is computing with T-symmetry. change 'true' to 'false' if you want to compute the whole BZ 
; true -> half of BZ, false -> whole BZ 
(define-param T-symmetry? true)

;choose one method to use: GGR or tetrahedron (Tr) method.
;the grid of two methods is different and both of two need to cover the complete BZ volumn (or half of BZ)
;the 'false' is GGR method. change 'false' to 'true' if you want to use Tr method.
;false -> GGR, true -> Tr
(define-param GGR-Tr? false)

;num-k is the inter number of k points along one dimension when sampling the BZ
(define-param num-k 10); num-k should be even number to make the mesh avoid Gamma point when using GGR method or T-symmetry reduction

;---------------------------------------------------
; computing setting -- accuracy and number of bands
;---------------------------------------------------

(set-param! resolution 16)
(set-param! mesh-size 2)
(set-param! num-bands 20)

;-----------------------------------User changing part ends-----------------------------------------

;the BZ ranges in [-kx,+kx]^3
(define kx (/ 1 2))

;k-shift-GGR is half of the distance between two k ajacent points in one dimension 
(define k-shift-GGR (/ kx (+ 2 num-k)))
(define k-shift 0)
(if GGR-Tr?
	(set! k-shift 0); Tr k points
	(set! k-shift k-shift-GGR); GGR k points
)

;map the k points in Brillouin zone
(define (flatten x)
    (cond ((null? x) '())
          ((not (pair? x)) (list x))
          (else (append (flatten (car x))
                        (flatten (cdr x))))))

(define kpts
(flatten
(map
(lambda (z)
        (map
        (lambda (y)
                (map
                (lambda (x) (vector3 x y z))
		(if T-symmetry?
			(interpolate (- (/ num-k 2) 1) (list (- k-shift 0) (- kx k-shift))); half of the BZ
			(interpolate num-k (list (- k-shift kx) (- kx k-shift))); the whole BZ 
		)
                )
       )
        (interpolate num-k (list (- k-shift kx) (- kx k-shift)))
        )
)
	(interpolate num-k (list (- k-shift kx) (- kx k-shift)));
)
))

;compute the frequency band in line or map
(define-param Zone? true)
(if Zone?
        (set! k-points (map (lambda (x) x) kpts))
	(set! k-points (interpolate k-int bandseq))
)

;run the program
(if GGR-Tr?
	(run); Tr method doesn't need group velocity
	(run display-group-velocities); GGR method needs group velocity
)
