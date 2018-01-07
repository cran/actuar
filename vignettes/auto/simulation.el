(TeX-add-style-hook
 "simulation"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "x11names")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("natbib" "round") ("babel" "english") ("helvet" "scaled=0.9") ("mathpazo" "sc") ("enumitem" "shortlabels") ("Sweave" "noae" "inconsolata")))
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
    "amsthm"
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
    "exampleautorefname"
    "FrameCommand")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:rmixture"
    "eq:mixture"
    "sec:rcompound"
    "eq:definition-S"
    "ex:comppois"
    "sec:simul"
    "eq:basic_model"
    "sec:simul:description"
    "sec:simul:usage")
   (LaTeX-add-environments
    '("exercice" LaTeX-env-args ["argument"] 0)
    '("emphbox" LaTeX-env-args ["argument"] 1)
    '("texample" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-bibliographies
    "actuar")
   (LaTeX-add-amsthm-newtheorems
    "example")
   (LaTeX-add-xcolor-definecolors
    "link"
    "url"
    "citation"
    "codebg"))
 :latex)

