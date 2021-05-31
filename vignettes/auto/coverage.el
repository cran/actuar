(TeX-add-style-hook
 "coverage"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names" "english")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("helvet" "scaled=0.9") ("mathpazo" "sc") ("enumitem" "shortlabels") ("Sweave" "noae")))
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
    "natbib"
    "babel"
    "helvet"
    "mathpazo"
    "enumitem"
    "Sweave"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    '("code" 1)
    '("pkg" 1)
    "D")
   (LaTeX-add-bibliographies
    "actuar")
   (LaTeX-add-xcolor-definecolors
    "link"
    "url"
    "citation"
    "codebg"))
 :latex)

