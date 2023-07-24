(TeX-add-style-hook
 "agg"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("report" "a4paper" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "report"
    "rep12"
    "times"
    "amsmath"
    "graphicx"
    "cite")
   (TeX-add-symbols)
   (LaTeX-add-labels
    "cha:theory"
    "sec:theory:sections")
   (LaTeX-add-bibliographies
    "abbrev"
    "maths"))
 :latex)

