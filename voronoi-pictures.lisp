(defpackage :voronoi-pictures
  (:use :cl)
  (:use :alexandria)
  (:use :opticl)
  (:use :lparallel)
  (:use :voronoi-pictures.app-utils)
  (:export :-main))

(in-package :voronoi-pictures)


(declaim (optimize (speed 3) (safety 3) (debug 3) (space 3) (compilation-speed 0)))

(defparameter *num-additions* 5)
(declaim (fixnum *num-additions*))

(defun range (min max)
  (declare (fixnum min max))
  (loop for i from min below max collecting i))

(defstruct v
  (x 0 :type number)
  (y 0 :type number)
  (points (make-set nil))
  (sum-color (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0))
             :type (simple-array fixnum (3)))
  (average-color (make-array 3 :element-type '(unsigned-byte 8) :initial-contents '(0 0 0))
                 :type (simple-array (unsigned-byte 8) (3)))
  (sse 0.0 :type number)
  (invalid t))

(defun copy-table (table)
  (let ((new-table (make-hash-table
                    :test (hash-table-test table)
                    :size (hash-table-size table))))
    (maphash #'(lambda(key value)
                 (setf (gethash key new-table) value))
             table)
    new-table))

(defun copy-voro (voro)
  (let* ((dimensions (array-dimensions voro))
         (new-array (make-array dimensions
                                :element-type (array-element-type voro)
                                :adjustable (adjustable-array-p voro)
                                :fill-pointer (and (array-has-fill-pointer-p voro)
                                                   (fill-pointer voro)))))
    (dotimes (i (length voro))
      (setf (row-major-aref new-array i)
            (copy-vo (row-major-aref voro i))))
    new-array))

(defun copy-vo (v)
  (make-v :x (v-x v)
          :y (v-y v)
          :points (copy-set (v-points v))
          :sum-color (copy-array (v-sum-color v))
          :average-color (copy-array (v-average-color v))
          :sse (v-sse v)
          :invalid (v-invalid v)))

(defun copy-set (set)
  (copy-table set))

(defun make-picture-array (img)
  (with-image-bounds (height width)
      img
    (let ((arr (make-array (list height width) :element-type 'fixnum )))
      (loop for i below height do
           (loop for j below width do
                (setf (aref arr i j) 0)))
      arr)))

(defun in-bound-point (x min max)
  (if (< x min)
      nil
      (if (> x max)
          nil
          x)))

(defmacro while (test &body body)
  `(loop until (not ,test)
      do (progn ,@body)))

(defun get-random-point-in-image (x-min x-max y-min y-max img)
  (declare (fixnum x-min x-max y-min y-max))
  (with-image-bounds (imh imw) img
    (let ((x -1)
          (y -1))
      (declare (fixnum x y))
      (while (not (and (in-bound-point x 0 imw)
                       (in-bound-point y 0 imh))) 
        (setf x (+ x-min (random (- x-max x-min))))
        (setf y (+ y-min (random (- y-max y-min)))))
      (list x y))))

(defun get-random-initial-point (img)
  (with-image-bounds (imh imw) img
    (get-random-point-in-image 0 imw 0 imh img)))

(defun initialize-voronoi-points (img &optional (num-additions *num-additions* ))
  (let ((v-arr (make-array 0 :element-type 'v :adjustable t :fill-pointer t)))
    (loop for i from 0 below num-additions do
         (destructuring-bind (x y)
             (get-random-initial-point img)
           (let ((v (make-v :x x :y y)))
             (vector-push-extend v v-arr))))
    v-arr))

(defun min-index (fn arr)
  (declare (function fn))
  (let ((mini 0)
        (minv (funcall fn (aref arr 0))))
    (loop for i from 0 below (length arr) do
         (let ((v (funcall fn (aref arr i))))
           (when (< v minv)
             (setf minv v)
             (setf mini i))))
    (values mini (aref arr mini))))

(defun square (x)
  (declare (number x))
  (the number (* x x)))


(defun dist-sq-v (x y v)
  (declare (fixnum x y))
  (+ (square (- (v-x v) x))
     (square (- (v-y v) y))))

(defun nearest-voronoi (x y oldmini v-arr)
  (declare (fixnum x y oldmini) ((vector v) v-arr))
  (let ((l (length v-arr))
        (mindist (dist-sq-v x y (aref v-arr oldmini)))
        (minind oldmini))
    (loop for i in (range (the fixnum (- l (the fixnum *num-additions*))) l) do 
         (let* ((v (aref v-arr i))
                (d  (dist-sq-v x y v)))
           (when (< d mindist)
             (setf mindist d)
             (setf minind i))))
    minind))

(defun sum-colors (sum r g b)
  (declare ((simple-array fixnum (3)) sum) ((unsigned-byte 8) r g b))
  (incf (aref sum 0) r)
  (incf (aref sum 1) g)
  (incf (aref sum 2) b))

(defun average-color (v)
  (declare (type v v))
  (let ((len (hash-table-count (v-points v))))
    (loop for i below (length (v-sum-color v)) do
         (setf (aref (v-average-color v) i) (round (aref (v-sum-color v) i) (if (= 0 len) 1 len))))))

(defun voronoi-bucket (arr v-arr minx maxx miny maxy &optional (first-run nil))
  (declare ((vector v) v-arr)
           (fixnum minx maxx miny maxy)
           ((simple-array fixnum (* *)) arr))
  (map nil (lambda (i)
             (declare (fixnum i))
             (loop for j from minx below maxx 
                do 
                  (let* ((oldmini (aref arr i j))
                         (oldv (aref v-arr oldmini))
                         (mini (nearest-voronoi j i oldmini v-arr))
                         (newv (aref v-arr mini)))
                    (when (or first-run (not (= mini oldmini)))
                      (when (not (v-invalid oldv))
                        (setf (v-invalid oldv) t)
                        (setf (v-sum-color oldv)
                              (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0))))
                      (when (not (v-invalid newv))
                        (setf (v-invalid newv) t)
                        (setf (v-sum-color newv)
                              (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0))))
                      (remhash (cons i j) (v-points oldv))
                      (set-key (cons i j) (v-points newv))
                      (setf (aref arr i j) mini)))))
       (range miny maxy)))

(defun voronoi-stat-collect (v-arr img)
  (declare (type 8-bit-rgb-image img) ((vector v) v-arr))
  (pmap nil
        (lambda (v)
          (when (v-invalid v)
            (loop for point being the hash-keys of (v-points v) do
                 (destructuring-bind (i . j) point
                   (multiple-value-bind (r g b) (pixel img i j)
                     (sum-colors (v-sum-color v) r g b))))))
        v-arr))

(defun image-convert (img)
  (with-image-bounds (height width) img
    (let ((new-img (make-8-bit-rgb-image height width)))
      (typecase img
        (8-bit-rgba-image
         (locally
             (declare (type 8-bit-rgba-image img))
           (loop for i below height do
                (loop for j below width do
                     (multiple-value-bind (r g b a) (pixel img i j)
                       (declare (ignore a))
                       (setf (pixel new-img i j) (values r g b)))))))
        (8-bit-rgb-image
         (locally
             (declare (type 8-bit-rgb-image img))
           (loop for i below height do
                (loop for j below width do
                     (multiple-value-bind (r g b) (pixel img i j)
                       (setf (pixel new-img i j) (values r g b))))))) )
      new-img)))

(defun fix-averages (voro)
  (loop for v across voro do
       (when (v-invalid v)
         (average-color v))))

(defun make-picture (voro-map voro img)
  (with-image-bounds (height width) img
    (let ((image (make-8-bit-rgb-image height width)))
      (loop for i below height do
           (loop for j below width do
                (let* ((v (aref voro (aref voro-map i j)))
                       (ac (v-average-color v)))
                  (setf (pixel image i j)
                        (values (aref ac 0)
                                (aref ac 1)
                                (aref ac 2))))))
      image)))

(defun fsquare (x)
  (declare (number x))
  (* x x))

(defmacro inc-err-channel (v refval val)
  `(incf (v-sse ,v) (fsquare (the number (- (the number ,refval)
                                            (the number ,val))))))


(defun calc-error (v ref-arr img)
  (declare (8-bit-rgb-image img))
  (let* ((vr (aref (v-average-color v) 0))
         (vg (aref (v-average-color v) 1))
         (vb (aref (v-average-color v) 2)))
    (multiple-value-bind (vl va vb) (rgb-lab vr vg vb)
      (loop for point being the hash-keys of (v-points v) do
           (multiple-value-bind (r g b) (pixel img (car point) (cdr point))
             (multiple-value-bind (l a b) (rgb-lab r g b)
               (inc-err-channel v vl l)
               (inc-err-channel v va a)
               (inc-err-channel v vb b)))))
    (let ((val (hash-table-count (v-points v))))
      (if (< val (square *num-additions*))
          (setf (v-sse v) 0.0)
          (setf (v-sse v) (/ (v-sse v) (sqrt val)))))))

(defun shift-average (v amount)
  (when (> amount 0.001)
    (let* ((acol (v-average-color v)))
      (multiple-value-bind (l a b) (rgb-lab (aref acol 0)
                                            (aref acol 1)
                                            (aref acol 2))
        
        (let ((l (+ l (- (* 100.0 amount) (random (* 200.0 amount)))))
              (a (+ a (- (* 100.0 amount) (random (* 200.0 amount)))))
              (b (+ b (- (* 100.0 amount) (random (* 200.0 amount))))))
          (multiple-value-bind (r g b) (lab-rgb l a b)
            (setf (aref acol 0) r)
            (setf (aref acol 1) g)
            (setf (aref acol 2) b)))))))

(defun shift-averages (voros amount)
  (pmap nil (lambda (v) (shift-average v amount)) voros))

(defun calc-errors (voro img)
  (pmap nil (lambda (v)
              (when (v-invalid v)
                (calc-error v nil img)
                (setf (v-invalid v) nil)))
        voro))

(defun find-worst-cell (voro)
  (min-index (lambda (v) (- (v-sse v))) voro))

(defun hash-table-to-array (table)
  (let ((vec (make-array (hash-table-count table))))
    (loop
       for key being the hash-keys of table 
       for i = 0 then (1+ i) do
         (setf (aref vec i) key))
    vec))

(defun split-cell (voro v &optional (num-additions *num-additions*))
  (let ((s (shuffle (hash-table-to-array (v-points v)))))
    (loop for i below (min num-additions (length s)) do
         (let ((p (aref s i)))
           (when (not (and (= (v-x v) (car p))
                           (= (v-y v) (cdr p))))
             (vector-push-extend (make-v :x (cdr p)
                                         :y (car p)) voro))))))

(defun split-lowest-cell (voro)
  (let ((v (nth-value 1 (find-worst-cell voro))))
    (split-cell voro v)
    v))

(defun minimum (seq) 
  (the fixnum (reduce #'min seq :initial-value (elt seq 0))))

(defun maximum (seq)
  (the fixnum (reduce #'max seq :initial-value (elt seq 0))))

(defun get-v-bounds (v)
  (let* ((points (hash-table-to-array (v-points v)))
         (xs (map 'list #'car points))
         (ys (map 'list #'cdr points))
         (mpi (min-index (lambda (p)
                           (the fixnum
                                (- ;; Invert the order, so get arg-max, not arg-min
                                 (+
                                  (square (the (signed-byte 32)
                                               (- (v-x v) (the fixnum (cdr p)))))
                                  (square (the (signed-byte 32)
                                               (- (v-y v) (the fixnum (car p)))))))))
                         points))
         (mp (aref points mpi))
         (d (+ 1
               (ceiling (isqrt (the (signed-byte 32)
                                    (+
                                     (square (the (signed-byte 32)
                                                  (- (v-x v) (the fixnum (cdr mp)))))
                                     (square (the (signed-byte 32)
                                                  (- (v-y v) (the fixnum (car mp))))))))))))
    (declare (fixnum d))
    (values (- (minimum xs) d) (+ (maximum xs) d)
            (- (minimum ys) d) (+ (maximum ys) d))))


(defun optimize-loop (voro img max color-shift &optional video)
  (let ((voro voro)
        (arr (make-picture-array img)))
    (with-image-bounds (height width) img
      (let ((minx 0) (maxx width)
            (miny 0) (maxy height))
        (voronoi-bucket arr voro minx maxx miny maxy t)
        (loop for i below max do
             (progn
               (format t "~a~a%        " #\return (/ (floor (* 10000 (/ i (* max 1.0)))) 100.0))
               (finish-output)
               (voronoi-bucket arr voro minx maxx miny maxy)
               (voronoi-stat-collect voro img)
               (fix-averages voro)                                
               (calc-errors voro img)
               (when video
                 (let ((pic (make-picture arr voro img)))
                   
                   (write-image-file (format nil "~a/video-~12,'0d.png" video i) pic)))
               (multiple-value-bind (nminx nmaxx nminy nmaxy)
                   (get-v-bounds (split-lowest-cell voro))
                 (setf minx (max 0 nminx))
                 (setf maxx (min width nmaxx))
                 (setf miny (max 0 nminy))
                 (setf maxy (min height nmaxy))))))
      (voronoi-bucket arr voro 0 width 0 height t)
      (voronoi-stat-collect voro img)
      (fix-averages voro)
      (when color-shift (shift-averages voro color-shift))
      (format t "~a100%        ~%" #\return)
      (values arr voro))))

(defmacro rgb-xyz-correct-channel (chan)
  `(setf ,chan (* 100 (if (> ,chan 0.04045)
                          (expt (/ (+ ,chan 0.055) 1.055) 2.4)
                          (/ ,chan 12.92)))))

(defmacro linear-combine (wr wg wb)
  `(+ (* vr ,wr) (* vg ,wg) (* vb ,wb)))

(defparameter *ref_x*  95.047)
(defparameter *inv_ref_x* (/ 1 *ref_x*))
(declaim (type number *inv_ref_x*))

(defparameter *ref_y* 100.000)
(defparameter *inv_ref_y* (/ 1 *ref_y*))
(declaim (type number *inv_ref_y*))

(defparameter *ref_z* 108.883)
(defparameter *inv_ref_z* (/ 1 *ref_z*))
(declaim (type number *inv_ref_z*))

(defparameter *inv255* (/ 1 255.0))
(declaim (type number *inv255*))

(defun rgb-xyz (r g b) 
  (declare (unsigned-byte r g b))
  (let ((vr (* r *inv255*))
        (vg (* g *inv255*))
        (vb (* b *inv255*)))
    (declare (number vr vg vb))
    (rgb-xyz-correct-channel vr)
    (rgb-xyz-correct-channel vg)
    (rgb-xyz-correct-channel vb)
    (values
     (linear-combine 0.4124 0.3576 0.1805)
     (linear-combine 0.2126 0.7152 0.0722)
     (linear-combine 0.0193 0.1192 0.9505))))

(defun clip (v mi ma) 
  (if (> v ma)
      ma
      (if (< v mi)
          mi
          v)))

(defmacro xyz-rgb-correct-channel (chan)
  `(progn (setf ,chan
                (round
                 (the number
                      (clip
                       (the number
                            (* 255.0 (if (> ,chan 0.0031308)
                                         (- (* 1.055 (expt ,chan (/ 2.4)))
                                            0.055)
                                         (* 12.92 ,chan))))
                       0.0
                       255.0))))))

(defun xyz-rgb (x y z)
  (declare (number x y z))
  (let ((vr (* x 0.01))
        (vg (* y 0.01))
        (vb (* z 0.01)))
    (declare (number vr vg vb))
    (let ((vr (linear-combine 3.2406 -1.5372 -0.4986))
          (vg (linear-combine -0.9692 1.8758 0.0415))
          (vb (linear-combine 0.0557 -0.2040 1.0570)))
      (xyz-rgb-correct-channel vr)
      (xyz-rgb-correct-channel vg)
      (xyz-rgb-correct-channel vb)
      (values vr vg vb))))


(defmacro xyz-lab-correct-channel (chan)
  `(setf ,chan (if (> ,chan 0.008856)
                   (expt ,chan (/ 1 3))
                   (+ (/ 16.0 116) (* ,chan 7.787)))))

(defun xyz-lab (x y z)
  (declare (number x y z))
  (let ((vx (* x *inv_ref_x*))
        (vy (* y *inv_ref_y*))
        (vz (* z *inv_ref_z*)))
    (xyz-lab-correct-channel vx)
    (xyz-lab-correct-channel vy)
    (xyz-lab-correct-channel vz)
    (values (- (* 116 vy) 16)
            (* 500 (- vx vy))
            (* 200 (- vy vz)))))

(defun cube (x)
  (* x (* x x)))

(defun lab-xyz (l a b)
  (let*  ((fy (/ (+ l 16.0) 116.0))
          (fz (- fy (/ b 200.0)))
          (fx (+ fy (/ a 500.0)))
          (fx3 (cube fx))
          (fz3 (cube fz)))
    (let ((k 903.3)
          (e 0.008856))
      (let ((x (if (> fx3 e)
                   fx3
                   (/ (- (* 116.0 fx) 16) k)))
            (y (if (> l (* k e))
                   (cube fy)
                   (/ l k)))
            (z (if (> fz3 e)
                   fz3
                   (/ (- (* 116.0 fz) 16) k))))
        (values
         (* x *ref_x*)
         (* y *ref_y*)
         (* z *ref_z*))))))


(defun rgb-lab (r g b)
  (multiple-value-bind (x y z) (rgb-xyz r g b)
    (xyz-lab x y z)))

(defun lab-rgb (l a b)
  (multiple-value-bind (x y z) (lab-xyz l a b)
    (xyz-rgb x y z)))

(defun voro-to-kd-points (voro) 
  (loop for i below (length voro) collecting
       (let ((v (aref voro i)))
         (make-array '(3)
                     :initial-contents (list (v-x v) (v-y v) i) :adjustable nil
                     :element-type 'number :fill-pointer nil))))

(defun sort-kd-voro (voro-kd)
  (let* ((x-a (copy-list voro-kd))
         (y-a (copy-list voro-kd))
         (x-a (sort x-a #'< :key (lambda (v) (aref v 0))))
         (y-a (sort y-a #'< :key (lambda (v) (aref v 1)))))
    (values x-a y-a)))

(defun make-set (list)
  (let ((set (make-hash-table :test 'equal :size 128)))
    (loop for x in list do
         (setf (gethash x set) t))
    set))

(defun key-in-set (key set)
  (gethash key set))

(defun set-key (key set)
  (setf (gethash key set) t))

(defun collect-axis (voro-set axis-sort)
  (let (res)
    (loop for point in axis-sort do
         (when (key-in-set point voro-set)
           (push point res)))
    (nreverse res)))

(defun split-list-at-point (list n)
  (declare (fixnum n))
  (if (= n 0)
      (values nil list)
      (loop repeat n
         for tail = list then (cdr tail)
         collecting (car tail) into head
         finally (return (values head (cdr tail))))))

(defun split-axis (axis-sort)
  (declare (list axis-sort))
  (multiple-value-bind (left right)
      (split-list-at-point axis-sort (floor (length axis-sort) 2))
    (values left (car right) (cdr right))))

(defpun make-kd-tree-helper (voro curr-ax other-ax)
  (when voro
    (multiple-value-bind (left med right) (split-axis voro)
      (plet ((left-set (make-set left)) 
             (right-set (make-set right)))
        (plet ((left (collect-axis left-set other-ax))
               (right (collect-axis right-set other-ax)))
          (plet ((lt (make-kd-tree-helper left  left curr-ax))
                 (rt (make-kd-tree-helper right right curr-ax)))
            (list med
                  lt
                  rt)))))))

(defun make-kd-tree-from-voro (voro)
  (let ((points (voro-to-kd-points voro)))
    (make-kd-tree points)))

(defun make-kd-tree (kd-points)
  (multiple-value-bind (x-ax y-ax) (sort-kd-voro kd-points)
    (make-kd-tree-helper x-ax x-ax y-ax)))

(defun dist (x y point)
  (declare (number x y))
  (sqrt (+ (square (- x (aref point 0)))
           (square (- y (aref point 1))))))

(defun nearest-neighbor-helper (kd-tree x y axis nearest nearest-dist)
  (declare (number x y))
  (when kd-tree
    (let* ((p (first kd-tree))
           (d (dist x y p))
           (val (* 1.0 (aref p axis))))
      (declare 
       (number d)
       (number val))
      (when (< d (the number (first nearest-dist)))
        (setf (first nearest) p)
        (setf (first nearest-dist) d))
      (if (= axis 0)
          (if (< x val)
              (progn (nearest-neighbor-helper (second kd-tree) x y 1 nearest nearest-dist)
                     (when (> (+ (the number (first nearest-dist)) x) val)
                       (nearest-neighbor-helper (third kd-tree) x y 1 nearest nearest-dist)))
              (progn (nearest-neighbor-helper (third kd-tree) x y 1 nearest nearest-dist)
                     (when (< (- x (the number (first nearest-dist))) val)
                       (nearest-neighbor-helper (second kd-tree) x y 1 nearest nearest-dist))))
          (if (< y val)
              (progn (nearest-neighbor-helper (second kd-tree) x y 0 nearest nearest-dist)
                     (when (> (+ (the number (first nearest-dist)) y) val)
                       (nearest-neighbor-helper (third kd-tree) x y 0 nearest nearest-dist)))
              (progn (nearest-neighbor-helper (third kd-tree) x y 0 nearest nearest-dist)
                     (when (< (- y (the number (first nearest-dist))) val)
                       (nearest-neighbor-helper (second kd-tree) x y 0 nearest nearest-dist))))))))

(defun nearest-neighbor (kd-tree x y)
  (declare (number x y))
  (let ((nearest (list (first kd-tree)))
        (nearest-dist (list (dist x y (first kd-tree)))))
    (nearest-neighbor-helper kd-tree x y 0 nearest nearest-dist)
    (aref (first nearest) 2)))


(defun open-image (file)
  (image-convert (read-image-file file)))

(defun runner (in-file out-file iterations &optional color-shift)
  (progn
    (setf lparallel:*kernel* (lparallel:make-kernel 8))
    (let* ((img (open-image in-file))
           (voro (initialize-voronoi-points img)))
      (multiple-value-bind (ar nvoro) (optimize-loop voro img iterations color-shift)
        (let ((lab-img (make-lab-img img)))
          (print (eval-voro nvoro img lab-img)))
        (let ((out-img (make-picture ar nvoro img)))
          (write-image-file out-file out-img))))))

(defun -main-old (&optional args)
  (if (not (>= (length args) 4))
      (format t "Usage: ~a <input-file> <output-file> <number-of-iterations> [color-shift] ~%" (and args (first args)))
      (runner (pathname (second args))
              (pathname  (third args))
              (parse-integer (fourth args))
              (if (fifth args)
                  (/ (parse-integer (fifth args)) 100.0)
                  0.0))))


(defun voronoi-diagram (v-arr minx maxx miny maxy)
  (declare ((vector v) v-arr)
           (fixnum minx maxx miny maxy))
  (let ((kd-tree (make-kd-tree-from-voro v-arr)))
    (map nil (lambda (i)
               (declare (fixnum i))
               (loop for j from minx below maxx 
                  do 
                    (let* ((mini (nearest-neighbor kd-tree j i))
                           (newv (aref v-arr mini)))
                      (when (not (v-invalid newv))
                        (setf (v-invalid newv) t)
                        (setf (v-sum-color newv)
                              (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0)))
                        (setf (v-points newv) (make-set (list))))
                      (set-key (cons i j) (v-points newv))
                      ;;(setf (aref arr i j) mini)
                      )))
         (range miny maxy))))

(defun voronoi-arr (arr v-arr minx maxx miny maxy)
  (declare ((vector v) v-arr)
           (fixnum minx maxx miny maxy))
  (let ((kd-tree (make-kd-tree-from-voro v-arr)))
    (map nil (lambda (i)
               (declare (fixnum i))
               (loop for j from minx below maxx 
                  do 
                    (let* ((mini (nearest-neighbor kd-tree j i)))
                      (setf (aref arr i j) mini)
                      )))
         (range miny maxy))))

(defun eval-perf ()
  (let* ((img (open-image "../../../test.png"))
         (lab-img (make-lab-img img))
         (voros (time (loop for i below 20 collecting
                           (initialize-voronoi-points img 1000)))))
    (time (pmap 'list
                (lambda (voro)
                  (eval-voro voro img lab-img))
                voros))))

(defun calc-sum-error (v lab-img)
  (declare ((simple-array number (* * 3)) lab-img))
  (let* ((vr (aref (v-average-color v) 0))
         (vg (aref (v-average-color v) 1))
         (vb (aref (v-average-color v) 2)))
    (multiple-value-bind (vl va vb) (rgb-lab vr vg vb)
      (loop for point being the hash-keys of (v-points v) do
           (let ((y (car point))
                 (x (cdr point)))
             (let ((l (aref lab-img y x 0))
                   (a (aref lab-img y x 1))
                   (b (aref lab-img y x 2)))
               (inc-err-channel v vl l)
               (inc-err-channel v va a)
               (inc-err-channel v vb b)))))
    (v-sse v)))

(defun make-lab-img (img)
  (declare (8-bit-rgb-image img))
  (with-image-bounds (height width) img
    (let ((lab-img (make-array (list height width 3) :element-type 'number)))
      (loop for i below height do
           (loop for j below width do
                (multiple-value-bind (r g b) (pixel img i j)
                  (multiple-value-bind (l a b) (rgb-lab r g b)
                    (setf (aref lab-img i j 0) l)
                    (setf (aref lab-img i j 1) a)
                    (setf (aref lab-img i j 2) b)))))
      lab-img)))

(defun calc-sum-errors (voro lab-img)
  (reduce #'+
          (pmap 'list (lambda (v)
                        (calc-sum-error v lab-img))
                voro)
          :initial-value 0))

(defun eval-voro (voro img lab-img)
  (loop for v across voro do
       (progn
         (setf (v-invalid v) nil)
         (setf (v-sum-color v)
               (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0)))
         (setf (v-sse v) 0.0)))
  (with-image-bounds (height width) img
    (voronoi-diagram voro 0 width 0 height))
  (voronoi-stat-collect voro img)
  (fix-averages voro)                                
  (calc-sum-errors voro lab-img))



(defun random? (&optional (prob 0.5))
  "Tosses a coin of prob probability of coming up heads,
then returns t if it's heads, else nil."
  (< (random 1.0) prob))

(defun generate-list (num function &optional no-duplicates)
  "Generates a list of size NUM, with each element created by
  (funcall FUNCTION).  If no-duplicates is t, then no duplicates
are permitted (FUNCTION is repeatedly called until a unique
new slot is created).  EQUALP is the test used for duplicates."
  (let (bag)
    (while (< (length bag) num) bag
           (let ((candidate (funcall function)))
             (unless (and no-duplicates
                          (member candidate bag :test #'equalp))
               (push candidate bag))))
    bag))

(defun random-pair (firsts seconds)
  (random-elt
   (mapcar #'cons
		   firsts
		   seconds)))

(defparameter *tournament-size* 5)
(defun tournament-select-one (population fitnesses)
  "Does one tournament selection and returns the selected individual."
  (let ((participants
		 (generate-list *tournament-size*
						(lambda ()
						  (random-pair population fitnesses)))))
    (copy-voro (first (reduce
                       (lambda (a b) 
                         (if (> (cdr a) (cdr b))
                             a
                             b))
                       participants)))))


(defun tournament-selector (num population fitnesses)
  "Does NUM tournament selections, and puts them all in a list"
  (loop repeat num
	 collect (tournament-select-one population fitnesses)))

(defun best-by (evaluator collection)
  (reduce
   (lambda (a b)
     (if (> (funcall evaluator a) (funcall evaluator b))
		 a
		 b))
   collection))

(defun best-pair (values scores)
  (let ((pairs (mapcar #'cons
					   values scores)))
	(best-by #'cdr pairs)))

(defun evolve (generations pop-size
               &key setup creator selector modifier evaluator printer writer)
  "Evolves for some number of GENERATIONS, creating a population of size
POP-SIZE, using various functions"
  ;; The functions passed in are as follows:
  ;;(SETUP)                     called at the beginning of evolution, to set up
  ;;                            global variables as necessary
  ;;(CREATOR)                   creates a random individual
  ;;(SELECTOR num pop fitneses) given a population and a list of corresponding fitnesses,
  ;;                            selects and returns NUM individuals as a list.
  ;;                            An individual may appear more than once in the list.
  ;;(MODIFIER ind1 ind2)        modifies individuals ind1 and ind2 by crossing them
  ;;                            over and mutating them.  Returns the two children
  ;;                            as a list: (child1 child2).  Nondestructive to
  ;;                            ind1 and ind2.
  ;;(PRINTER pop fitnesses)     prints the best individual in the population, plus
  ;;                            its fitness, and any other interesting statistics
  ;;                            you think interesting for that generation.
  ;;(EVALUATOR individual)      evaluates an individual, and returns its fitness.
  ;;Pop will be guaranteed to be a multiple of 2 in size.
  ;;
  ;; HIGHER FITNESSES ARE BETTER
  
  ;; your function should call PRINTER each generation, and also print out or the
  ;; best individual discovered over the whole run at the end, plus its fitness
  ;; and any other statistics you think might be nifty.
  
    ;;; IMPLEMENT ME
  ;; setup pre-first-gen
  (funcall setup)
  (let* ((current-population
		  (loop repeat pop-size
			 collecting (funcall creator)))
		 (current-fitnesses
		  (pmap 'list evaluator current-population))
		 (best '(() . ()))
		 (gen 0)
		 greatest-generation)
    ;; now do stuff for each generation
    (loop repeat generations
       do
		 (progn
           (when (= (mod gen 5) 0)
             (sb-ext:gc :full t))
		   (incf gen)
		   (let ((moms (funcall selector (/ pop-size 2) ; choose which individuals reproduce
								current-population
								current-fitnesses))
				 (dads (funcall selector (/ pop-size 2)
								current-population
								current-fitnesses)))
			 (setf current-population ; make children and cull
				   (mapcan modifier moms dads)))
		   (setf current-fitnesses ; recalculate fitnesses
				 (pmap 'list evaluator current-population))
           (let ((new-best
                  (best-pair current-population current-fitnesses)))
             (if (or (not (cdr best))
                     (> (cdr new-best) (cdr best)))
                 (progn
                   (setf best new-best)
                   (setf greatest-generation gen)
                   (format t "Current Generation: ~a~%" gen)
                   (format t "  Current Best Fitness: ~a~%" (cdr best))
                   (funcall writer (car best) gen)

                   )))))
    ;; write up the entire run
	(format t "~%~%Best Result:~%")
    (format t "Greatest Generation: ~D~%" greatest-generation)
    (funcall printer (list (car best)) (list (cdr best)))
    (car best)))

;; img lab-img
(defun gen-runner ()
  (let* ((img (open-image "../pres/leaf0.png"))
         (lab-img (make-lab-img img)))
    (let ((voro (evolve 20 25 :setup (lambda ())
                        :creator (lambda () (initialize-voronoi-points img 5000))
                        :selector #'tournament-selector
                        :modifier (lambda (mom dad) (modifier mom dad img))
                        :evaluator (lambda (voro) (- (eval-voro voro img lab-img)))
                        :printer (lambda (_ fits) (print (first fits)))
                        :writer (lambda (voro gen)
                                  (when nil
                                    (with-image-bounds (height width) img
                                      (let ((arr (make-picture-array img)))
                                        (voronoi-arr arr voro 0 width 0 height)
                                        (let ((pic (make-picture arr voro img)))
                                          (write-image-file (format nil "voroout/~12,'0d.png" gen)
                                                            pic)))))))))
      voro)))


(defun voro-modifier (vo width height)
  (if (random? 0.05)
      (progn
        (setf (v-x vo) (random width))
        (setf (v-y vo) (random height)))
      (progn
        (setf (v-x vo) (+ (* 3 (gaussian-random))
                          (v-x vo)))
        (setf (v-y vo) (+ (* 3 (gaussian-random))
                          (v-y vo)))))
  vo)

(defun vl-modifier (voro img)
  (if (random 0.02)
      (initialize-voronoi-points img 5000)
      (progn
        (when (random? 0.5)
          (let ((i (random (length voro))))
            (swap (aref voro i) (aref voro (1- (length voro))))
            (vector-pop voro)
            (split-cell voro (aref voro (find-worst-cell voro)) 1)))
        (with-image-bounds (height width) img
          (map 'vector (lambda (vo) (voro-modifier vo width height)) voro)))))

(defun modifier (mom dad img)
  (let ((c1 (copy-voro mom))
        (c2 (copy-voro dad)))
    (loop for i below (length mom) do
         (when (random? 0.05)
           (setf (aref c1 i) (aref dad i))
           (setf (aref c2 i) (aref mom i))))
    (vl-modifier c1 img)

    (vl-modifier c2 img)
    (list c1 c2)))

(defun -main (&optional args)
  (setf *random-state* (make-random-state t))
  (setf lparallel:*kernel* (lparallel:make-kernel 8))
  (let* ((voro (gen-runner))
         (img (open-image "../pres/leaf0.png"))
         (arr (make-picture-array img)))
    (with-image-bounds (height width) img
      (voronoi-arr arr voro 0 width 0 height))
    (let ((pic (make-picture arr voro img)))
      (write-image-file "out.png" pic))))


(defun find-best-cell (voro)
  (min-index (lambda (v) (v-sse v)) voro))

(defmacro swap (a1 a2)
  (let ((tmp (gensym)))
    `(let ((,tmp ,a1))
       (setf ,a1 ,a2)
       (setf ,a2 ,tmp))))
