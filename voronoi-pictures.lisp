(defpackage :voronoi-pictures
  (:use :cl)
  (:use :opticl)
  (:use :voronoi-pictures.app-utils)
  (:export :-main))

(in-package :voronoi-pictures)


(declaim (optimize (speed 3) (safety 0) (debug 0)))

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


(defstruct v
  (x 0 :type fixnum)
  (y 0 :type fixnum)
  (points (make-array 0 :element-type 'list :adjustable t :fill-pointer t) :type (vector list))
  (sum-color (make-array 3 :element-type 'fixnum :initial-contents '(0 0 0))
             :type (simple-array fixnum (3)))
  (average-color (make-array 3 :element-type '(unsigned-byte 8) :initial-contents '(0 0 0))
                 :type (simple-array (unsigned-byte 8) (3)))
  (sse 0.0 :type single-float))

(defun make-picture-array (img)
  (with-image-bounds (height width)
      img
    (make-array (list height width) :element-type 'fixnum)))

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

(defun nearest-voronoi (x y v-arr)
  (declare (fixnum x y))
  (min-index (lambda (v) (declare (type v v))
                     (the fixnum (+ (square (- x (v-x v))) (square (- y (v-y v))))))
             v-arr))

(defun sum-colors (sum col)
  (declare (vector sum))
  (loop for i below (length sum) do
       (incf (aref sum i) (elt col i))))

(defun average-color (v)
  (declare (type v v))
  (let ((len (length (v-points v))))
    (loop for i below (length (v-sum-color v)) do
         (setf (aref (v-average-color v) i) (floor (aref (v-sum-color v) i) (if (= 0 len) 1 len))))))

(defun voronoi-bucket (v-arr img)
  (declare (type 8-bit-rgb-image img) ((vector v) v-arr))
  (let ((arr (make-picture-array img)))
    (with-image-bounds (height width) img
      (loop for i below height
         do (loop for j below width 
               do 
                 (let ((mini (nearest-voronoi j i v-arr)))
                   (vector-push-extend (list i j) (v-points (aref v-arr mini)))
                   (sum-colors (v-sum-color (aref v-arr mini)) (multiple-value-list (pixel img i j)))
                   (setf (aref arr i j) mini)))))
    arr))



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
       (average-color v)))

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
  (loop for point across (v-points v) do
       (multiple-value-bind (r g b) (pixel img (first point) (second point))
         (inc-err-channel v 0 r)
         (inc-err-channel v 1 g)
         (inc-err-channel v 2 b)))
  (let ((val (length (v-points v))))
    (if (= 0 val)
        (setf (v-sse v) 0.0)
        (setf (v-sse v) (/ (v-sse v) (sqrt val))))))

(defun calc-errors (voro img)
  (loop for v across voro do
       (calc-error v img)))

(defun find-worst-cell (voro)
  (min-index (lambda (v) (- (v-sse v))) voro))

(defun split-cell (voro v)
  (loop for i below *num-additions* do
       (let ((point (aref (v-points v) (random (length (v-points v))))))
         (vector-push-extend (make-v :x (second point)
                                     :y (first point)) voro))))
(defun split-lowest-cell (voro)
  (let ((v (nth-value 1 (find-worst-cell voro))))
    (split-cell voro v)))

(defun optimize-voros (voro)
  (let ((new-arr (make-array 0 :element-type 'v :adjustable t :fill-pointer t)))
    (loop for i below (length voro) do
         (let ((v (aref voro i)))
           (when (> (length (v-points v)) 0)
             (vector-push-extend (make-v :x (v-x v)
                                         :y (v-y v)
                                         :sse (v-sse v)
                                         :average-color (v-average-color v)
                                         :sum-color (v-sum-color v)
                                         :points (v-points v)) new-arr))))
    new-arr))

(defun reset-voros (voro)
  (declare ((vector v) voro))
  (loop for i below (length voro) do
       (let ((v (aref voro i)))
         (setf (aref voro i) (make-v :x (v-x v)
                                     :y (v-y v))))))

(defparameter *num-additions* 5)

(defun optimize-loop (voro img)
  (let ((voro voro))
    (loop for i below 250 do
         (progn
           (format t "~a~%" i)           
           (reset-voros voro)
           (voronoi-bucket voro img)
           (setf voro (optimize-voros voro))
           (fix-averages voro)
           (calc-errors voro img)
           (split-lowest-cell voro)))
    (let ((ar (voronoi-bucket voro img)))
      (fix-averages voro)
      (values ar voro))))

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

(defun open-image ()
  (image-convert (read-png-file "/home/jsvlrt/Pictures/up-throw.png")))


(defun -main (&optional args)
  (format t "~a~%" "I don't do much yet"))


