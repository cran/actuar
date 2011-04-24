(TeX-add-style-hook "simulation"
 (lambda ()
    (LaTeX-add-environments
     "example")
    (LaTeX-add-labels
     "sec:introduction")
    (TeX-add-symbols
     '("code" 1)
     '("pkg" 1)
     '("proglang" 1))
    (TeX-run-style-hooks
     "Sweave"
     "noae"
     "paralist"
     "mathpazo"
     "sc"
     "helvet"
     "babel"
     "english"
     "natbib"
     "round"
     "hyperref"
     "color"
     "amsthm"
     "amsmath"
     "inputenc"
     "utf8"
     "fontenc"
     "T1"
     "latex2e"
     "art10"
     "article")))

