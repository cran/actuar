(TeX-add-style-hook
 "risk"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names" "english")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("babel" "english") ("helvet" "scaled=0.9") ("mathpazo" "sc") ("enumitem" "shortlabels") ("Sweave" "noae" "inconsolata")))
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
    "framed"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    '("code" 1)
    '("pkg" 1)
    '("mat" 1)
    '("VAR" 1)
    '("E" 1)
    "VaR"
    "CTE"
    "FrameCommand")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:collective-risk-model"
    "eq:definition-S"
    "eq:cdf-S"
    "eq:convolution-formula"
    "sec:discretization"
    "eq:discretization:upper"
    "eq:discretization:lower"
    "eq:discretization:midpoint"
    "eq:discretization:unbiased"
    "fig:discretization-methods"
    "sec:aggregate"
    "eq:normal-approximation"
    "eq:np2-approximation"
    "fig:Fs"
    "eq:VaR"
    "eq:CTE"
    "fig:Fs-comparison"
    "sec:ruin-model"
    "eq:definition-surplus"
    "eq:definition-ruin"
    "eq:definition-S(t)"
    "eq:definition-wait"
    "sec:adjustment-coefficient"
    "eq:definition-adjcoef"
    "eq:definition-adjcoef-ind"
    "eq:definition-adjcoef-prop"
    "eq:definition-adjcoef-xl"
    "fig:adjcoef"
    "sec:ruin"
    "eq:ruin-cramer-lundberg"
    "eq:prob-ruin:cramer-lundberg"
    "eq:prob-ruin:sparre:pi+"
    "eq:eq:prob-ruin:sparre:Q"
    "eq:kronecker-sum"
    "fig:prob-ruin"
    "sec:beekman"
    "eq:beekman:prob-ruin"
    "eq:beekman:p"
    "eq:beekman:H"
    "eq:beekman:prob-ruin-long"
    "fig:beekman:prob-ruin")
   (LaTeX-add-bibliographies
    "actuar")
   (LaTeX-add-xcolor-definecolors
    "link"
    "url"
    "citation"
    "codebg"))
 :latex)

