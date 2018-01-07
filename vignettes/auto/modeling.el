(TeX-add-style-hook
 "modeling"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("babel" "english") ("numprint" "autolanguage") ("helvet" "scaled=0.9") ("mathpazo" "sc") ("enumitem" "shortlabels") ("Sweave" "noae" "inconsolata")))
   (add-to-list 'LaTeX-verbatim-environments-local "Verbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "fontenc"
    "inputenc"
    "amsmath"
    "amsthm"
    "natbib"
    "babel"
    "numprint"
    "helvet"
    "mathpazo"
    "booktabs"
    "enumitem"
    "Sweave"
    "framed"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    '("code" 1)
    '("pkg" 1)
    '("mat" 1)
    '("VAR" 1)
    '("E" 1)
    "LAS"
    "exampleautorefname"
    "FrameCommand")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:grouped-data"
    "ex:grouped.data-1"
    "ex:grouped.data-2"
    "fig:histogram"
    "eq:ecdf"
    "ex:ogive"
    "fig:ogive"
    "sec:data-sets"
    "sec:empirical-moments"
    "fig:elev"
    "sec:minimum-distance"
    "ex:mde"
    "sec:coverage"
    "tab:coverage"
    "ex:coverage"
    "eq:pdf-YP")
   (LaTeX-add-environments
    '("exercice" LaTeX-env-args ["argument"] 0)
    '("emphbox" LaTeX-env-args ["argument"] 1)
    '("texample" LaTeX-env-args ["argument"] 0)
    "example"
    "rem")
   (LaTeX-add-bibliographies
    "actuar"))
 :latex)

