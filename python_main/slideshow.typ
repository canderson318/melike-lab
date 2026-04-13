// find outputs/ -type f \( -name "*.pdf" -o -name "*.png" \) | sed "s|$(pwd)||" | sed 's/.*/"&",/' >> slideshow.typ 
// typst compile slide-show.typ 

#let slide(body) = {
  pagebreak(weak: true)
  align(center + horizon)[#body]
}

#set page(margin: 0.5in, width: 8.5in, height:11in)


#let files = (
"outputs/03/SM002/param_against_interval_start.pdf",
"outputs/03/SM002/param_corr_heamap.pdf",
"outputs/03/SM002/param_pairplot.pdf",
"outputs/04/param_against_HROD_boxplots.pdf",
"outputs/04/param_against_deltaday_boxplots.pdf"
)


#for ( f ) in files{
    slide()[
         #stack(
            align(center)[#image(f)],
            v(-1em),
            align(bottom+right, text(size: .5em, strong[#f])),
        )
    ]
}
