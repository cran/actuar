(TeX-add-style-hook
 "credibility"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("babel" "english") ("helvet" "scaled=0.9") ("mathpazo" "sc") ("Sweave" "noae" "inconsolata")))
   (add-to-list 'LaTeX-verbatim-environments-local "Verbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
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
    "Sweave"
    "framed"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    '("code" 1)
    '("pkg" 1)
    '("E" 1)
    "pt"
    "FrameCommand")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:hachemeister"
    "sec:hierarchical"
    "eq:hierarchical:premiums"
    "eq:s2"
    "eq:Xbzzw"
    "eq:ac-BG"
    "eq:bc-BG"
    "eq:ac-Ohl"
    "eq:bc-Ohl"
    "eq:at"
    "eq:bt"
    "eq:Xzzw"
    "sec:buhlmann"
    "eq:a-hat"
    "eq:a-tilde"
    "sec:regression"
    "fig:state4"
    "fig:state4:2")
   (LaTeX-add-environments
    '("exercice" LaTeX-env-args ["argument"] 0)
    '("emphbox" LaTeX-env-args ["argument"] 1)
    '("texample" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-bibliographies
    "actuar")
   (LaTeX-add-xcolor-definecolors
    "link"
    "url"
    "citation"
    "codebg"))
 :latex)

