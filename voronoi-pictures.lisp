(defpackage :voronoi-pictures
  (:use :cl)
  (:use :alexandria)
  (:use :opticl)
  (:use :lparallel)
  (:use :voronoi-pictures.app-utils)
  (:export :-main))

(in-package :voronoi-pictures)


(declaim (optimize (speed 3) (safety 0) (debug 0) (space 0) (compilation-speed 0)))

(defun do-image-invert ()
  (let ((img (read-png-file "test.png")))
    (typecase img
      (8-bit-rgba-image
       (locally
           (declare (type 8-bit-rgba-image img))
         (with-image-bounds (height width)
             img
           (time
            (loop for i below height
               do (loop for j below width 
                     do 
                       (multiple-value-bind (r g b a)
                           (pixel img i j)
                         (declare (type (unsigned-byte 8) r g b a))
                         (setf (pixel img i j)
                               (values (- 255 r) g b a))))))))))
    (write-png-file "testinv.png" img)))

(defun range (min max)
  (declare (fixnum min max))
  (loop for i from min below max collecting i))

(defstruct v
  (x 0 :type fixnum)
  (y 0 :type fixnum)
  (points (make-set nil))
  (sum-color (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0))
             :type (simple-array fixnum (3)))
  (average-color (make-array 3 :element-type '(unsigned-byte 8) :initial-contents '(0 0 0))
                 :type (simple-array (unsigned-byte 8) (3)))
  (sse 0.0 :type single-float)
  (invalid t))

(defun make-picture-array (img)
  (with-image-bounds (height width)
      img
    (make-array (list height width) :element-type 'fixnum )))

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

(defun initialize-voronoi-points (img)
  (let ((v-arr (make-array 0 :element-type 'v :adjustable t :fill-pointer t)))
    (loop for i from 0 below *num-additions* do
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
  (declare (fixnum x))
  (the fixnum (* x x)))

(defun nearest-voronoi (x y kd-tree)
  (declare (fixnum x y))
  (aref (the (simple-array fixnum (3)) (nearest-neighbor kd-tree x y)) 2))

(defun sum-colors (sum r g b)
  (declare ((simple-array fixnum (3)) sum) ((unsigned-byte 8) r g b))
  (incf (aref sum 0) r)
  (incf (aref sum 1) g)
  (incf (aref sum 2) b))

(defun average-color (v)
  (declare (type v v))
  (let ((len (hash-table-count (v-points v))))
    (loop for i below (length (v-sum-color v)) do
         (setf (aref (v-average-color v) i) (floor (aref (v-sum-color v) i) (if (= 0 len) 1 len))))))

(defun voronoi-bucket (arr v-arr minx maxx miny maxy &optional (first-run nil))
  (declare ((vector v) v-arr)
           ((simple-array fixnum (* *)) arr))
  (let ((kd-tree (make-kd-tree v-arr)))
    (map nil (lambda (i)
               (loop for j from minx below maxx 
                  do 
                    (let* ((mini (nearest-voronoi j i kd-tree))
                           (newv (aref v-arr mini))
                           (oldmini (aref arr i j))
                           (oldv (aref v-arr oldmini)))
                      (when (not (v-invalid oldv))
                        (setf (v-invalid oldv) t)
                        (setf (v-sum-color oldv) (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0))))
                      (when (not (v-invalid newv))
                        (setf (v-invalid newv) t)
                        (setf (v-sum-color newv) (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0))))
                      (when (not (or first-run (= mini oldmini)))
                        (remhash (list i j) (v-points oldv))
                        (set-key (list i j) (v-points newv))
                        (setf (aref arr i j) mini)))))
         (range miny maxy))))


(defun voronoi-stat-collect (v-arr img)
  (declare (type 8-bit-rgb-image img) ((vector v) v-arr))
  (pmap nil
        (lambda (v)
          (when (v-invalid v)
            (loop for point being the hash-keys of (v-points v) do
                 (destructuring-bind (i j) point
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

(defmacro inc-err-channel (v chan var)
  `(incf (v-sse ,v) (square (- (aref (v-average-color ,v) ,chan) ,var))))


(defun calc-error (v img)
  (declare (8-bit-rgb-image img))
  (loop for point being the hash-keys of (v-points v) do
       (multiple-value-bind (r g b) (pixel img (first point) (second point))
         (inc-err-channel v 0 r)
         (inc-err-channel v 1 g)
         (inc-err-channel v 2 b)))
  (let ((val (hash-table-count (v-points v))))
    (if (< val (square *num-additions*))
        (setf (v-sse v) 0.0)
        (setf (v-sse v) (/ (v-sse v) (sqrt val))))))

(defun calc-errors (voro img)
  (pmap nil (lambda (v)
              (when (v-invalid v)
                (calc-error v img)
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

(defun split-cell (voro v)
  (let ((s (shuffle (hash-table-to-array (v-points v)))))
    (loop for i below (min *num-additions* (length s)) do
         (let ((p (aref s i)))
           (when (not (and (= (v-x v) (second p))
                           (= (v-y v) (first p))))
             (vector-push-extend (make-v :x (second p)
                                         :y (first p)) voro))))))

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
         (xs (map 'list #'second points))
         (ys (map 'list #'first points))
         (mpi (min-index (lambda (p) (the fixnum
                                          (- (+ (square (the (signed-byte 32) (- (v-x v) (the fixnum (second p)))))
                                                (square (the (signed-byte 32) (- (v-y v) (the fixnum (first p))))))))) points))
         (mp (aref points mpi))
         (d (ceiling (isqrt (the (signed-byte 32)
                                 (+ (square (the (signed-byte 32) (- (v-x v) (the fixnum (second mp)))))
                                    (square (the (signed-byte 32) (- (v-y v) (the fixnum (first mp)))))))) 2)))
    (declare (fixnum d))
    (values (- (minimum xs) d) (+ (maximum xs) d)
            (- (minimum ys) d) (+ (maximum ys) d))))

(defparameter *num-additions* 5)

(defun clear-voro (voro)
  (loop for i below (length voro) do
       (setf (aref voro i) (make-v :x (v-x (aref voro i))
                                   :y (v-y (aref voro i))))))

(defun optimize-loop (voro img max)
  (let ((voro voro)
        (arr (make-picture-array img)))
    (with-image-bounds (height width) img
      (let ((minx 0) (maxx width)
            (miny 0) (maxy height))
        (voronoi-bucket arr voro minx maxx miny maxy t)
        (clear-voro voro)
        (loop for i below max do
             (progn
               (format t "~a~a%        " #\return (/ (floor (* 10000 (/ i (* max 1.0)))) 100.0))
               (finish-output)
               (voronoi-bucket arr voro minx maxx miny maxy)
               (voronoi-stat-collect voro img)
               (fix-averages voro)
               (calc-errors voro img)
               (multiple-value-bind (nminx nmaxx nminy nmaxy)
                   (get-v-bounds (split-lowest-cell voro))
                 (setf minx (max 0 nminx))
                 (setf maxx (min width nmaxx))
                 (setf miny (max 0 nminy))
                 (setf maxy (min height nmaxy))))))
      (voronoi-bucket arr voro 0 width 0 height)
      (voronoi-stat-collect voro img)
      (fix-averages voro)
      (format t "~a100%        ~%" #\return)
      (values arr voro))))

(defun rgb-hsv (r g b)
  (let* ((min (* (/ 255.0) (min r g b)))
         (max (* (/ 255.0) (max r g b)))
         (v max)
         (delta (- max min))
         s h)
    (when (< delta 0.01)
      (setf delta 0.01))
    (if (< max 0.00001)
        (progn (setf s 0)
               (setf h -1)
               (return-from rgb-hsv (values 0 0 0)))
        (progn
          (setf s (/ delta max))
          (if (= r max)
              (setf h (/ (- g b) delta))
              (if (= g max)
                  (setf h (+ 2 (/ (- b r) delta)))
                  (setf h (+ 4 (/ (- r g) delta)))))
          (setf h (* h 60))))
    (if (< h 0)
        (incf h 360.0))
    (values (coerce (max 0 (min 255 (floor (* (/ 255.0 360.0) h)))) '(unsigned-byte 8))
            (coerce (max 0 (min 255 (floor (* 255.0 s)))) '(unsigned-byte 8))
            (coerce (max 0 (min 255 (floor (* 255.0 v)))) '(unsigned-byte 8)))))

(defun hsv-rgb (h s v)
  (let (i f p q tt
          (h (* (/ 360.0 255.0) h))
          (s (* (/ 255.0) s))
          (v (* (/ 255.0) v)))
    (setf h (/ h 60.0))
    (setf i (floor h))
    (setf f (- h i))
    (setf p (* v (- 1 s)))
    (setf q (* v (- 1 (* s f))))
    (setf tt (* v (- 1 (* s (- 1 f)))))
    (multiple-value-bind (r g b)
        (case i
          (0 (values v tt p))
          (1 (values q v p))
          (2 (values p v tt))
          (3 (values p q v))
          (4 (values tt p v))
          (5 (values v p q)))
      (values (floor (* r 255.0))
              (floor (* g 255.0))
              (floor (* b 255.0))))))

(defun convert-image-to-hsv (img)
  (with-image-bounds (height width) img
    (loop for i below height do
         (loop for j below width do
              (multiple-value-bind (r g b) (pixel img i j)
                (setf (pixel img i j) (rgb-hsv r g b))))))
  img)

(defun convert-image-to-rgb (img)
  (with-image-bounds (height width) img
    (loop for i below height do
         (loop for j below width do
              (multiple-value-bind (h s v) (pixel img i j)
                (setf (pixel img i j) (hsv-rgb h s v))))))
  img)




(defun voro-to-kd-points (voro) 
  (loop for i below (length voro) collecting
       (let ((v (aref voro i)))
         (make-array '(3)
                     :initial-contents (list (v-x v) (v-y v) i) :adjustable nil
                     :element-type 'fixnum :fill-pointer nil))))

(defun sort-kd-voro (voro-kd)
  (let* ((x-a (copy-list voro-kd))
         (y-a (copy-list voro-kd))
         (x-a (sort x-a #'< :key (lambda (v) (aref v 0))))
         (y-a (sort y-a #'< :key (lambda (v) (aref v 1)))))
    (values x-a y-a)))

(defun make-set (list)
  (let ((set (make-hash-table :test 'equal)))
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

(defun make-kd-tree (voro)
  (let ((points (voro-to-kd-points voro)))
    (multiple-value-bind (x-ax y-ax) (sort-kd-voro points)
      (make-kd-tree-helper x-ax x-ax y-ax))))

(defun dist (x y point)
  (declare (fixnum x y) ((simple-array fixnum (3)) point))
  (the single-float (sqrt (+ (the (signed-byte 32) (square (the (signed-byte 32) (- x (aref point 0)))))
                             (the (signed-byte 32) (square (the (signed-byte 32) (- y (aref point 1)))))))))

(defun nearest-neighbor-helper (kd-tree x y axis nearest nearest-dist)
  (declare (fixnum x y))
  (when kd-tree
    (let* ((p (first kd-tree))
           (d (dist x y p))
           (val (* 1.0 (aref p axis))))
      (declare ((simple-array fixnum (3)) p)
               (single-float d)
               (single-float val))
      (when (< d (the single-float (first nearest-dist)))
        (setf (first nearest) p)
        (setf (first nearest-dist) d))
      (if (= axis 0)
          (if (< x val)
              (progn (nearest-neighbor-helper (second kd-tree) x y 1 nearest nearest-dist)
                     (when (> (+ (the single-float (first nearest-dist)) x) val)
                       (nearest-neighbor-helper (third kd-tree) x y 1 nearest nearest-dist)))
              (progn (nearest-neighbor-helper (third kd-tree) x y 1 nearest nearest-dist)
                     (when (< (- x (the single-float (first nearest-dist))) val)
                       (nearest-neighbor-helper (second kd-tree) x y 1 nearest nearest-dist))))
          (if (< y val)
              (progn (nearest-neighbor-helper (second kd-tree) x y 0 nearest nearest-dist)
                     (when (> (+ (the single-float (first nearest-dist)) y) val)
                       (nearest-neighbor-helper (third kd-tree) x y 0 nearest nearest-dist)))
              (progn (nearest-neighbor-helper (third kd-tree) x y 0 nearest nearest-dist)
                     (when (< (- y (the single-float (first nearest-dist))) val)
                       (nearest-neighbor-helper (second kd-tree) x y 0 nearest nearest-dist))))))))

(defun nearest-neighbor (kd-tree x y)
  (declare (fixnum x y))
  (let ((nearest (list (first kd-tree)))
        (nearest-dist (list (dist x y (first kd-tree)))))
    (nearest-neighbor-helper kd-tree x y 0 nearest nearest-dist)
    (the (simple-array fixnum (3)) (first nearest))))

(defun open-image (file)
  (image-convert (read-image-file file)))

(defun runner (in-file out-file iterations)
  (progn
    (setf lparallel:*kernel* (lparallel:make-kernel 8))
    (let* ((img (open-image in-file))
           (voro (initialize-voronoi-points img)))
      (multiple-value-bind (ar nvoro) (optimize-loop voro img iterations)
        (let ((out-img (make-picture ar nvoro img)))
          (write-image-file out-file out-img))))))

(defun -main (&optional args)
  (if (not (= (length args) 4))
      (format t "Usage: ~a <input-file> <output-file> <number-of-iterations> ~%" (and args (first args)))
      (runner (pathname (second args))
              (pathname  (third args))
              (parse-integer (fourth args)))))
