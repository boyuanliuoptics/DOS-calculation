;-----------------------------------User changing part begins-----------------------------------------

;-------------------
; structure setting
;-------------------

; A triangular lattice of air rods in dielectric 13.

;set the primitive lattice
(set! geometry-lattice (make lattice (size 1 1 no-size)
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))

;define a couple of parameters about the structure
(set! geometry (list (make cylinder
                       (center 0 0 0) (radius 0.48) (height infinity)
                       (material (make dielectric (epsilon 1))))))

(set! default-material (make dielectric (epsilon 13)))
;-----------------------------------------------------
; computing setting -- for drawing band stucture line
;-----------------------------------------------------

;switch for computing the band in the line 
;if you want only want to calculate DOS, change 'true' to any other word. It will be detected in the shell file.
;switch-> true

;high symmetry points in the band
(define G (vector3 0 0 0))
(define M (vector3 0 0.5 0))
(define K (vector3 (/ -3) (/ 3) 0))

;sequence of band structure in line
(define bandseq (list G M K G))

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
(define-param num-k 30); num-k should be even number to make the mesh avoid Gamma point when using GGR method or T-symmetry reducement

;---------------------------------------------------
; computing setting -- accuracy and number of bands
;---------------------------------------------------

(set-param! resolution 16)
(set-param! mesh-size 2)
(set-param! num-bands 8)

;-----------------------------------User changing part ends-----------------------------------------

;the BZ ranges in [-kx,+kx]^2
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
	(list 0)
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
	(run-tm); Tr method doesn't need group velocity
	(run-tm display-group-velocities); GGR method needs group velocity
)
