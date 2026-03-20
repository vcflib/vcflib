;; Generate a Singularity container with:
;;
;;   guix pack -f squashfs -S /bin=bin -m guix-container-manifest.scm
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
