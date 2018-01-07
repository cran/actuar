(TeX-add-style-hook
 "actuar"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("babel" "english") ("helvet" "scaled=0.9") ("mathpazo" "sc") ("Sweave" "noae" "inconsolata")))
   (add-to-list 'LaTeX-verbatim-environments-local "Verbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "fontenc"
    "inputenc"
    "amsmath"
    "natbib"
    "babel"
    "helvet"
    "mathpazo"
    "Sweave"
    "framed"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    '("pkg" 1)
    "FrameCommand")
   (LaTeX-add-labels
    "sec:introduction")
   (LaTeX-add-environments
    '("exercice" LaTeX-env-args ["argument"] 0)
    '("emphbox" LaTeX-env-args ["argument"] 1)
    '("texample" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-bibliographies))
 :latex)

