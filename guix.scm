;; To use this file to build HEAD of vcflib:
;;
;;   guix build -f guix.scm
;;   guix build -L . vcflib-local-htslib-git  # test local htslib build - fails because of path rewrites
;;   guix build -L . vcflib-static-git --tune=native
;;
;; To get a development container (emacs shell will work)
;;
;;   guix shell -C -D -F -f guix.scm
;;
;;   cmake  -DCMAKE_BUILD_TYPE=Debug -DOPENMP=OFF -DASAN=ON ..
;;   make -j 12
;;   ctest .
;;
;; debug example
;;
;;   env LD_LIBRARY_PATH=$GUIX_ENVIRONMENT/lib gdb --args vcfallelicprimitives -m ../samples/10158243.vcf
;;
;;   guix shell -L . -C -D -F vcflib-static-git
;;
;; support other (external) zig compiler
;;
;; To bring in a recent zig compiler I do something like
;;
;;   guix shell -C -D -f guix.scm --expose=/home/wrk/opt/zig-linux-x86_64-0.11.0-dev.987+a1d82352d/=/zig
;;
;; and add /zig to the PATH. E.g.
;;
;;   export PATH=/zig:$PATH


(define-module (guix)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix build-system cmake)
  #:use-module (guix download)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module (guix packages)
  #:use-module (guix utils)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages autotools)
  #:use-module (gnu packages base)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages build-tools)
  #:use-module (gnu packages check)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages curl)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages gdb)
  #:use-module (gnu packages haskell-xyz) ; pandoc for help files
  #:use-module (gnu packages llvm)
  #:use-module (gnu packages parallel)
  #:use-module (gnu packages perl)
  #:use-module (gnu packages perl6)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages python)
  #:use-module (gnu packages python-xyz) ; for pybind11
  #:use-module (gnu packages ruby)
  #:use-module (gnu packages time)
  #:use-module (gnu packages tls)
  #:use-module (gnu packages zig)
  #:use-module (ice-9 popen)
  #:use-module (ice-9 rdelim)
  #:use-module (srfi srfi-1)
  )

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define %version
  (read-string (open-pipe "git describe --always --tags --long|tr -d $'\n'" OPEN_READ)))

(define-public vcflib-git
  (package
    (name "vcflib-git")
    (version %version)
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    ;; (arguments
    ;; `(#:tests? #t
    ;;   #:configure-flags
    ;;   ,#~(list
    ;;       ;; "-DBUILD_OPTIMIZED=ON"       ;; we don't use the standard cmake optimizations
    ;;       "-DCMAKE_BUILD_TYPE=Debug"))) ;; to optimize use guix --tune=march-type (e.g. --tune=native)
    (inputs
     ;; ("gcc" ,gcc-13)       ;; test against latest - won't build python bindings
     (list
       fastahack    ;; dev version not in Debian
       htslib ;; disable to test local build, also disable tabixpp in that case
       pandoc ; for man pages
       perl
       python
       python-pytest
       pybind11
       ruby ; for man pages
       smithwaterman
       tabixpp
       time ; for tests
       wfa2-lib ; alternative:  cmake  -DCMAKE_BUILD_TYPE=Debug -DWFA_GITMODULE=ON -DZIG=ON ..
       xz
       zig))
    (native-inputs
     `(("pkg-config" ,pkg-config)))
    (home-page "https://github.com/vcflib/vcflib/")
    (synopsis "Library for parsing and manipulating VCF files")
    (description "Vcflib provides methods to manipulate and interpret
sequence variation as it can be described by VCF.  It is both an API for parsing
and operating on records of genomic variation as it can be described by the VCF
format, and a collection of command-line utilities for executing complex
manipulations on VCF files.")
    (license license:expat)))

(define-public vcflib-test-one-git
  "Build and test one target"
  (package
    (inherit vcflib-git)
    (name "vcflib-test-one-git")
    (arguments
     `(#:tests? #t
       #:configure-flags
       ,#~(list
           ;; "-DBUILD_OPTIMIZED=ON"       ;; we don't use the standard cmake optimizations
           "-DCMAKE_BUILD_TYPE=Debug") ;; to optimize use guix --tune=march-type (e.g. --tune=native)
       #:phases
       ,#~(modify-phases %standard-phases
                         (delete 'install)
                         (replace 'build
                                  (lambda* (#:key make-flags #:allow-other-keys)
                                    (invoke "make" "-j 12" "pyvcflib")))
                         (replace 'check
                                  (lambda* (#:key tests? #:allow-other-keys)
                                    (when tests?
                                      (invoke "ctest" "-R" "pyvcflib" "--verbose")))))))))


(define-public vcflib-local-htslib-git
  "Test embedded htslib - part of tabixpp submodule"
  (package
    (inherit vcflib-git)
    (name "vcflib-local-htslib-git")
    (arguments
     `(#:tests? #f ;; tests don't work when running build directly
       #:configure-flags
       ,#~(list
           ;; "-DBUILD_OPTIMIZED=ON"       ;; we don't use the standard cmake optimizations
           "-DCMAKE_BUILD_TYPE=Generic"))) ;; to optimize use guix --tune=march-type (e.g. --tune=native)
    (inputs
     (modify-inputs (package-inputs vcflib-git)
                    (delete "htslib" "tabixpp")
                    (prepend
                     autoconf  ;; htslib build requirements
                     automake
                     libdeflate
                     openssl
                     curl)))))


;; ==== The following is for static binary builds using gcc - used mostly for deployment ===

;; Guix does not come with a static version of libdeflate
(define-public libdeflate-static
  (package
    (inherit libdeflate)
    (name "libdeflate-static")
    (version "1.19")
    (arguments
     (list #:configure-flags
           #~(list "-DLIBDEFLATE_BUILD_STATIC_LIB=YES"
                   "-DLIBDEFLATE_BUILD_TESTS=YES")))))

;; A minimal static version of htslib that does not depend on curl and openssl. This
;; reduces the number of higher order dependencies in static linking.
(define-public htslib-static
  (package
    (inherit htslib)
    (name "htslib-static")
    (version "1.19")
    (source (origin
            (method url-fetch)
            (uri (string-append
                  "https://github.com/samtools/htslib/releases/download/"
                  version "/htslib-" version ".tar.bz2"))
            (sha256
             (base32
              "0dh79lwpspwwfbkmllrrhbk8nkvlfc5b5ib4d0xg5ld79w6c8lc7"))))
    (arguments
     (substitute-keyword-arguments (package-arguments htslib)
       ((#:configure-flags flags ''())
        ''())))
    (inputs
     (list bzip2 xz))))

(define-public vcflib-static-git
  "Optimized for latest AMD architecture build and static deployment.
These binaries can be copied to HPC."
  (package
    (inherit vcflib-git)
    (name "vcflib-static-git")
    (arguments
     `(#:tests? #f
       #:configure-flags
       ,#~(list
           "-DBUILD_STATIC=ON"
           ;; "-DZIG=OFF"
           ;; "-DBUILD_OPTIMIZED=ON"    ;; we don't use the standard cmake optimizations
           "-DCMAKE_BUILD_TYPE=Generic" ;; to optimize use guix --tune=march-type (e.g. --tune=native)
           "-DCMAKE_INSTALL_RPATH=")))   ; force cmake static build and do not rewrite RPATH
    (inputs
     (modify-inputs (package-inputs vcflib-git)
                    (prepend
                     `(,bzip2 "static")
                     `(,zlib "static")
                     `(,xz "static")
                     libdeflate-static
                     htslib-static)))))


vcflib-git
