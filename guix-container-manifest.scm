;; Create a Singularity container with:
;;
;;   guix pack -f squashfs -S /bin=bin -L . -m guix-container-manifest.scm   # Singularity
;;
;; Create a Docker container with:
;;
;;   guix pack -f docker   -S /bin=bin -L . -m guix-container-manifest.scm   # Docker                
;;
;; Run with produced container:
;;
;;   singularity exec vcflib.squashfs vcfwave --help
;;   singularity exec vcflib.squashfs vcfcreatemulti --help

(use-modules (guix))

(packages->manifest
 (list vcflib-git
       (@ (gnu packages bash) bash-minimal)
       (@ (gnu packages base) coreutils-minimal)))
