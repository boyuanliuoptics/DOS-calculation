(define-param b 1)
(define-param c 1)
(define pi (* 4 (atan 1))) ; 3.14159...
(set! geometry-lattice (make lattice
(basis-size (/ (sqrt 3) 2) (/ (sqrt 3) 2) (/ (sqrt 3) 2))
(basis1 -1 1 1)
(basis2 1 -1 1)
(basis3 1 1 -1)))
;b1=(0 1 1); b2=(1 0 1); b3=(1 1 0)
; Corners of the irreducible Brillouin zone for the bcc lattice,
(define G  (vector3 0 0 0))
(define H (vector3 0.5 -0.5 0.5))
(define P (vector3 0.75 -0.25 -0.25))
(define N (vector3 0.5 0. -0.5))
;----------------
(define (flatten x)
    (cond ((null? x) '())
          ((not (pair? x)) (list x))
          (else (append (flatten (car x))
                        (flatten (cdr x))))))
(define-param num-k 10) 
(define kx 1)
(define ky kx)
(define kz kx)

;kcenter is in the reciprocal units!
;the kcenter denotes an arbitrary shift of mesh in BZ to avoid von Hove singularity like Gamma
(define kcenter (cartesian->reciprocal (vector3 0.07 0.08 0.09)))

;k points sampling on map/line
; sampling is to cover all the first BZ using the subcells so there is a term of (- kx (/ 1 (+ 2 num-k)))
(define kpts
(map cartesian->reciprocal
(flatten
(map
(lambda (z)
        (map
        (lambda (y)
                (map
                (lambda (x) (reciprocal->cartesian (vector3 x y z)))
                (interpolate num-k (list 0 (- kx (/ 1 (+ 2 num-k)))))
                )
       )
        (interpolate num-k (list 0 (- ky (/ 1 (+ 2 num-k)))))
        )
)
        (interpolate num-k (list 0 (- kz (/ 1 (+ 2 num-k)))))
))
))

(define-param k-int 30)
(define-param Zone? false)
(if Zone?
        (set! k-points (map (lambda (x) (vector3+ kcenter x)) kpts))
        (set! k-points (interpolate k-int (list H G N P G)))
)
;----------------


; define a couple of parameters (which we can set from the command-line)
(define-param eps 16) 

(define-param epsxx eps) ; the dielectric
(define-param epsyy eps)
(define-param epszz eps)

(define epsbg 1) ; the dielectric
(define-param isoval1 -1.1) ;isovalue to choose packing fraction
(define-param isoval2 isoval1) ;isovalue to choose packing fraction;(define diel (make dielectric (epsilon epsc)))

; The Epsilon function that returns the dielectric constant as a function of position
; returns epsc if the value of the level-set function V > isoval, or else returns 1

(define FB
(lambda (h k l x y z)
(* (cos (* 2 pi (+ h k l) 0.25)) (+ 
	(* (sin (* 2 pi (+ (* h x) (/ l 4)))) (sin (* 2 pi (+ (* k y) (/ h 4)))) (sin (* 2 pi (+ (* l z) (/ k 4)))))
	(* (sin (* 2 pi (+ (* h y) (/ l 4)))) (sin (* 2 pi (+ (* k z) (/ h 4)))) (sin (* 2 pi (+ (* l x) (/ k 4)))))
	(* (sin (* 2 pi (+ (* h z) (/ l 4)))) (sin (* 2 pi (+ (* k x) (/ h 4)))) (sin (* 2 pi (+ (* l y) (/ k 4)))))
))))

(define (eps-func p-lattice)
(let* ((p (lattice->cartesian p-lattice))
(x (vector3-x p))
(y (vector3-y p))
(z (vector3-z p)))
(if 
(< (FB 1 1 0 x (/ y b) (/ z c)) isoval1)
	(make dielectric-anisotropic (epsilon-diag epsxx epsyy epszz) (epsilon-offdiag 0 0 0))
	(if (< (FB 1 1 0 (- x) (- (/ y b)) (- (/ z c))) isoval2)
		(make dielectric-anisotropic (epsilon-diag epsxx epsyy epszz) (epsilon-offdiag 0 0 0))
		(make dielectric (epsilon epsbg)))
)
))

(set! default-material (make material-function (material-func eps-func)))

(set-param! resolution 16)
(set-param! mesh-size 2)
(set-param! num-bands 10)

(run display-group-velocities)

