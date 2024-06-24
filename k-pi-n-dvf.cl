;; K(pi,n) DISCRETE VECTOR FIELD FOR PI A CYCLIC GROUP




;; Auxiliar function for the bar-basis iterative process


(defun isog-list (list degr isog)
               (if (eq (length list) 0)
                   nil
               (if (eq (length list) 1)
                   (cons (funcall isog degr (first list)) '())
                 (cons (funcall isog degr (first list)) (isog-list (rest list) degr isog)))))



;; K(pi,1)^c BASIS 

;; Finite case: Zp

(defun aux-k1zp-basis-cc (p degr)
               (if (eq degr 0)
                   nil
                 (if (eq degr 1)
                     '(1)
                   (if (> degr 1)
                       (if (evenp degr)
                           (cons (- p 1) (aux-k1zp-basis-cc p (- degr 1)))
                         (cons 1 (aux-k1zp-basis-cc p (- degr 1))))))))

(defun basis-k1zp-cc (p)
               (flet ((rslt (degr)
                            (if (>= degr 0)
                                (cons (aux-k1zp-basis-cc p degr) '())
                              nil)))
                 #'rslt))


;; Infinite case: Z

(defun basis-k1z-cc ()
               (flet ((rslt (degr)
                            (if (eq degr 0)
                                (cons nil '())
                              (if (eq degr 1)
                                  (cons '(1) '())
                                    nil))))
                 #'rslt))


;; K(pi,n)^c CRITICAL COMPLEX


;; The following functions provide the nth iteration of the vector field algorithm
;; in the cases Zp and Z.

;; Iterated DVF

(defun iterated-knzp-vf (p n)
                (declare (type fixnum p))
                (declare (type fixnum n))
                (if (< n 1)
                    (error "n must be greater than zero.")
                  (if (eq n 1)
                      (kzp1-vf p)
                    (if (eq n 2)
                        (let(
                              (k1zp (k-zp p 1))
                              (k1zp-vf (kzp1-vf p)))
                          (cs-vf-g k1zp k1zp-vf))
                      (if (> n 2)
                          (cs-vf-g (k-zp p (- n 1))(iterated-knzp-vf p (- n 1))))))))


(defun iterated-knz-vf (n)
                (declare (type fixnum n))
                (if (< n 1)
                    (error "n must be greater than zero.")
                  (if (eq n 1)
                      (kz-vf)
                    (if (eq n 2)
                        (let(
                              (k1z (k-z-1))
                              (k1z-vf (kz-vf)))
                          (cs-vf-g k1z k1z-vf))
                      (if (> n 2)
                          (cs-vf-g (k-z (- n 1)) (iterated-knz-vf (- n 1))))))))


;; Iterated critical basis

(defun basis-cc-knzp (p n)
                (declare (type fixnum p))
                (declare (type fixnum n))
                (flet ((rslt (degr)
                             (declare (type fixnum degr))
                             (if (< n 1)
                                 (error "n must be greater than zero.")
                               (if (eq n 1)
                                   (funcall (basis-k1zp-cc p) degr)
                                 (if (eq n 2)
                                     (let*(
                                           (k1zp (k-zp p 1))
                                           (isog-k2zp (g-cs-isog-gnrt k1zp)))
                                       (isog-list (funcall (bar-basis (basis-k1zp-cc p)) degr) degr isog-k2zp))
                                   (if (> n 2)
                                       (let*(
                                             (kn1zp (k-zp p (- n 1)))
                                             (isog-kn1zp (g-cs-isog-gnrt kn1zp)))
                                         (isog-list (funcall (bar-basis (basis-cc-knzp p (- n 1))) degr) degr isog-kn1zp))))))))
                  #'rslt))


(defun basis-cc-knz (n)
                (declare (type fixnum n))
                (flet ((rslt (degr)
                             (declare (type fixnum degr))
                             (if (< n 1)
                                 (error "n must be greater than zero.")
                               (if (eq n 1)
                                   (funcall (basis-k1z-cc) degr)
                                 (if (eq n 2)
                                     (let*(
                                           (k1z (k-z 1))
                                           (isog-k2z (g-cs-isog-gnrt k1z)))
                                       (isog-list (funcall (bar-basis (basis-k1z-cc)) degr) degr isog-k2z))
                                   (if (> n 2)
                                       (let*(
                                             (kn1z (k-z (- n 1)))
                                             (isog-kn1z (g-cs-isog-gnrt kn1z)))
                                         (isog-list (funcall (bar-basis (basis-cc-knz (- n 1))) degr) degr isog-kn1z))))))))
                  #'rslt))

;; critical complex function



(defun knzp-cc (p n)
               (declare (type fixnum p))
               (declare (type fixnum n))
               (if (< n 1)
                   (error "n must be greater than zero.")
                 (let(
                       (knzp (k-zp p n))
                       (knzp-vf (iterated-knzp-vf p n)))
                   (build-chcm
                    :cmpr (cmpr knzp)
                    :basis (basis-cc-knzp p n)
                    :bsgn (crtc-complex-dffr knzp knzp-vf)
                    :intr-dffr (crtc-complex-dffr knzp knzp-vf)
                    :strt :gnrt
                    :orgn `(critical complex of k(Z ,p ,n))))))

(defun knz-cc (n)
               (declare (type fixnum n))
               (if (< n 1)
                   (error "n must be greater than zero.")
                 (let(
                       (knz (k-z n))
                       (knz-vf (iterated-knz-vf n)))
                   (build-chcm
                    :cmpr (cmpr knz)
                    :basis (basis-cc-knz n)
                    :bsgn (first (funcall (basis-cc-knz n) 0))
                    :intr-dffr (crtc-complex-dffr knz knz-vf)
                    :strt :gnrt
                    :orgn `(critical complex k(Z ,n))))))
