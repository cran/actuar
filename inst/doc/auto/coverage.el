(TeX-add-style-hook "coverage"
 (lambda ()
    (TeX-add-symbols
     "D")
    (TeX-run-style-hooks
     "Sweave"
     "amsmath"
     "latex2e"
     "art10"
     "article")))

