//
// find outputs/ -type f \( -name "*.pdf" -o -name "*.png" \) | sed "s|$(pwd)||"  > slideshow_paths.txt
// | sed 's/.*/"&",/'
// typst compile slide-show.typ 
//

#let slide(body) = {
  pagebreak(weak: true)
  align(center + horizon)[#body]
}

#set page(margin: 0.5in, width: 8.5in, height:11in)


#let files = read("slideshow_paths.txt").split("\n").slice(0,-1)

#for ( f ) in files{
    slide()[
         #stack(
            align(center)[#image(f)],
            v(-1em),
            align(bottom+right, text(size: .5em, strong[#f])),
        )
    ]
}
