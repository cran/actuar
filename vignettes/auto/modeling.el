(TeX-add-style-hook
 "modeling"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names" "english")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("numprint" "autolanguage") ("helvet" "scaled=0.9") ("mathpazo" "sc") ("enumitem" "shortlabels") ("Sweave" "noae" "inconsolata")))
   (add-to-list 'LaTeX-verbatim-environments-local "Verbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
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
   (LaTeX-add-bibliographies
    "actuar")
   (LaTeX-add-amsthm-newtheorems
    "example"
    "rem")
   (LaTeX-add-xcolor-definecolors
    "link"
    "url"
    "citation"
    "codebg"))
 :latex)

