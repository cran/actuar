(TeX-add-style-hook
 "distributions"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names" "english")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("helvet" "scaled=0.92") ("mathpazo" "sc") ("enumitem" "shortlabels") ("Sweave" "noae" "inconsolata")))
   (add-to-list 'LaTeX-verbatim-environments-local "Verbatim")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
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
    "microtype"
    "helvet"
    "mathpazo"
    "booktabs"
    "enumitem"
    "Sweave"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    '("samp" 1)
    '("code" 1)
    '("pkg" 1)
    '("mat" 1)
    '("VAR" 1)
    '("E" 1)
    "LAS")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:continuous"
    "tab:continuous"
    "sec:phase-type"
    "eq:Markov-transition-matrix"
    "eq:cdf-phtype"
    "eq:matrix-exponential"
    "sec:discrete"
    "tab:discrete"
    "eq:mixture"
    "sec:pig"
    "eq:pig:px"
    "eq:bessel_k"
    "eq:pig:px:recursive"
    "sec:special-integrals"
    "eq:pgamma"
    "eq:pbeta"
    "eq:gammainc"
    "eq:gammainc:apos"
    "eq:expint"
    "eq:betaint"
    "eq:betaint:2"
    "sec:implementation"
    "sec:app:continuous"
    "sec:app:continuous:transformed-beta"
    "sec:app:continuous:transformed-gamma"
    "sec:app:continuous:other"
    "sec:app:phase-type"
    "sec:app:discrete"
    "sec:app:discrete:a-b-0"
    "sec:app:discrete:zt"
    "sec:app:discrete:zm"
    "sec:app:discrete:pig")
   (LaTeX-add-environments
    '("texample" LaTeX-env-args ["argument"] 0)
    '("exercice" LaTeX-env-args ["argument"] 0)
    '("emphbox" LaTeX-env-args ["argument"] 1)
    '("bibexample" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-bibliographies
    "actuar")
   (LaTeX-add-xcolor-definecolors
    "link"
    "url"
    "citation"
    "codebg"))
 :latex)

