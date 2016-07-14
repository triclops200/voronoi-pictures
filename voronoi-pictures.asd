;;;; voronoi-pictures.asd

(asdf:defsystem #:voronoi-pictures
    :description "INSERT PROJECT DESCRIPTION HERE"
    :author "INSERT PROJECT AUTHOR HERE"
    :license "Modified BSD License"
    :serial t
    :depends-on (:opticl)
    :pathname "./"
    :components ((:file "app-utils")
                 (:file "voronoi-pictures")
                 ))

