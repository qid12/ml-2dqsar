(TeX-add-style-hook
 "MultiQASR_AutomationPhDMeet_160515"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "mathserif" "fleqn")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("babel" "english") ("tdclock" "timeinterval=1")))
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "beamerthemeshadow"
    "pgf"
    "pgfarrows"
    "pgfnodes"
    "pgfautomata"
    "pgfheaps"
    "amsmath"
    "amssymb"
    "amsfonts"
    "color"
    "xcolor"
    "graphicx"
    "multimedia"
    "manfnt"
    "babel"
    "algorithmic"
    "algorithm"
    "setspace"
    "tdclock"
    "tikz"
    "times")
   (TeX-add-symbols
    "toprule"
    "midrule"
    "botrule")
   (LaTeX-add-amsthm-newtheorems
    "jiashe"))
 :latex)

