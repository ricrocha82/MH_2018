# https://rpubs.com/Koundy/71792

# add new fonts
#install.packages('showtext', dependencies = TRUE)
#suppressPackageStartupMessages(library(showtext))
#showtext::showtext_auto()

# https://fonts.google.com/
#font_add_google( "Roboto Condensed", "RobotoCondensed-Regular")
#font_add_google("IBM Plex Serif", "IBM Plex Serif")
# Check the current search path for fonts
#font_paths()  

# syntax: font_add(family = "<family_name>", regular = "/path/to/font/file")
# font_add("Comic Sans MS", "comic.ttf")
#font_families()

library(extrafont)
# import fonts - only once
# font_import()
# load fonts - every session
loadfonts(device = "win", quiet = TRUE)

# find the name of a font you need for the family parameter of element_text
#library(tidyverse)
#fonttable() %>% filter(if_any(everything(), ~str_detect(.,"Comic")))

#fonttable() %>% dplyr::filter(str_detect(fontfile,"IBM")) %>% dplyr::pull(FamilyName)
#fonttable() %>% dplyr::filter(str_detect(fontfile,"Roboto")) %>% dplyr::pull(FamilyName)

# Automatically use showtext for new devices
#showtext_auto()
# Turn off if no longer needed
#showtext_auto(FALSE)

library(hrbrthemes)
library(ggthemes)
# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Using ggthemr package https://github.com/Mikata-Project/ggthemr
# set ggthemr theme
#ggthemr("<theme name>") 
# plot your existing figure with the new theme
#plt
# to remove all ggthemr effects later:
#ggthemr_reset()


# colors
# scales::show_col(c("#AD6F3B", "indianred3","#673770", "#8569D5", "#5E738F","coral4",
# "lightskyblue4", "firebrick4", "rosybrown1",
# "#CBD588", "orange","#DA5724", "#508578","#652926", "#C84248", "#D1A33D", "#8A7C64", "#599861","#5F7FC7","#CD9BCD"))
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#paletteer::paletteer_c("viridis::viridis", n= 20)
#barplot(c(2,5,4,5,6), col=c("burlywood", "deepskyblue2",'brown1',"forestgreen","darkmagenta" ))

#---------------------

theme_Publication_1 <- function(base_size = 14,
                                strip_text_size = 11,
                                strip_text_margin = 5,
                                subtitle_size = 13,
                                subtitle_margin = 10,
                                plot_title_size = 16,
                                plot_title_margin = 10,
                                title_axis_size = rel(1), # relative to the default
                                axis_text_size = 11,
                                ...) {
  
  ret <- ggplot2::theme_minimal(base_family = "helvetica",
                                base_size = base_size, ...)
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    # face = 'bold',
    margin = ggplot2::margin(b = strip_text_margin),
    family = "helvetica"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "helvetica"
  )
  ret$plot.title <- ggplot2::element_text(
    face = "bold",
    hjust = 0.5, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "helvetica"
  )
  ret$axis.title <- ggplot2::element_text(
    face = "bold",
    size = title_axis_size,
    family = "helvetica"
  )
  ret$axis.text <- ggplot2::element_text(
    face = "bold", 
    size = axis_text_size,
    family = "helvetica" 
  )
  ret
  
}

theme_Publication <- function(base_size=14, base_family="helvetica",...) {
  library(grid)
  library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.border = element_rect(fill = "transparent",colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(face = "bold",angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
            axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
          #  plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold", size = 11)
    ))
 
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# using the theme_article function form egg package

theme_Publication_2 <- function(...){

  egg::theme_article() +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(face = "bold", size = 12),
          axis.title.x = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 11), 
          axis.text.y = element_text(face = "bold", size = 11),
          strip.text = element_text(face="bold", size = 11),
          strip.background = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.margin = unit(0, "cm"),
          legend.title = element_text(face="italic", size = 12),
          legend.text = element_text(size = 11),
          text = element_text(family = "Garamond", color = "grey20"),
          legend.background = element_blank(),
          plot.margin=unit(c(10,5,5,5),"mm"))
 
  
}



#' Minimal ggplot2 theme using the Roboto Condensed and Roboto Bold fonts
#'
#' @param base_size base font size
#' @param strip_text_size,strip_text_margin plot strip text size and margin
#' @param subtitle_size,subtitle_margin plot subtitle size and margin
#' @param plot_title_size,plot_title_margin plot title size and margin
#' @param ... Other arguments passed to \code{theme_minimal}
#'
#' @details The Roboto Condensed and Roboto Bold fonts are both Google fonts;
#' they can be found at \url{https://fonts.google.com/specimen/Roboto+Condensed}
#' and \url{https://fonts.google.com/specimen/Roboto}. These fonts must be
#' installed locally on your computer for this theme to work.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     labs(title = "A Lovely Plot",
#'          subtitle = "What can the subtitle tell us?") +
#'     theme_roboto()
#'}
#'
#' @export

# for windows
theme_roboto <- function(base_size = 11,
                         strip_text_size = 12,
                         strip_text_margin = 5,
                         subtitle_size = 13,
                         subtitle_margin = 10,
                         plot_title_size = 16,
                         plot_title_margin = 10,
                         ...) {
  # Automatically use showtext for new devices
  
  ret <- ggplot2::theme_minimal(base_family = "Roboto",
                                base_size = base_size, ...)
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    margin = ggplot2::margin(b = strip_text_margin),
    family = "Roboto",
    face = 'bold'
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "Roboto",
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "Roboto",
    face = 'bold'
  ) 
  ret
  # Turn off if no longer needed
 
}




#' Minimal ggplot2 theme using the IBM Plex Sans fonts
#'
#' @param base_size base font size
#' @param strip_text_size,strip_text_margin plot strip text size and margin
#' @param subtitle_size,subtitle_margin plot subtitle size and margin
#' @param plot_title_size,plot_title_margin plot title size and margin
#' @param ... Other arguments passed to \code{theme_minimal}
#'
#' @details The IBM Plex fonts are open source and can be found at
#' \url{https://ibm.github.io/type/}. These fonts must be installed locally on
#' your computer for this theme to work.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     labs(title = "A Lovely Plot",
#'          subtitle = "What can the subtitle tell us?") +
#'     theme_plex()
#'
#' ggplot(diamonds, aes(carat, price, color = clarity)) +
#'     geom_point(alpha = 0.7) +
#'     facet_wrap(~cut) +
#'     labs(title = "A Lovely Plot",
#'          subtitle = "What can the subtitle tell us?") +
#'          theme_plex()
#'
#'}
#'
#' @export
theme_plex <- function(base_size = 11,
                       strip_text_size = 12,
                       strip_text_margin = 5,
                       subtitle_size = 13,
                       subtitle_margin = 10,
                       plot_title_size = 16,
                       plot_title_margin = 10,
                       ...) {

  ret <- ggplot2::theme_minimal(base_family = "IBM Plex Sans",
                                base_size = base_size, ...)
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    margin = ggplot2::margin(b = strip_text_margin),
    family = "IBM Plex Sans"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "IBM Plex Sans Sans"
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "IBM Plex Sans Sans"
  ) 
  ret
 
}


theme_Publication_3 <- function(base_size = 12,
                       strip_text_size = 12,
                       strip_text_margin = 5,
                       subtitle_size = 13,
                       subtitle_margin = 10,
                       plot_title_size = 16,
                       plot_title_margin = 10,
                       ...) {

  ret <- ggplot2::theme_minimal(base_family = "Comic Sans MS",
                                base_size = base_size, ...) +
    theme(panel.border = element_rect(fill = "NA",colour = "grey80")
    )
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
   # face = 'bold',
    margin = ggplot2::margin(b = strip_text_margin),
    family = "Comic Sans MS"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "Comic Sans MS"
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "Comic Sans MS"
  ) 
  ret

}


theme_depth_profile <- function(base_size = 12,
                                strip_text_size = 12,
                                strip_text_margin = 5,
                                subtitle_size = 13,
                                subtitle_margin = 10,
                                plot_title_size = 16,
                                plot_title_margin = 10,
                                ...) {
  
  ret <- theme_clean(base_family = "Comic Sans MS",
                                base_size = base_size, ...) +
    theme(plot.background = element_blank(),
          legend.background = element_blank()
    )
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    # face = 'bold',
    margin = ggplot2::margin(b = strip_text_margin),
    family = "Comic Sans MS"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "Comic Sans MS"
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "Comic Sans MS"
  ) 
  ret
  
}


# SC theme

theme_sc <- function(base_size = 12,
                     strip_text_size = 12,
                     strip_text_margin = 5,
                     subtitle_size = 13,
                     subtitle_margin = 10,
                     plot_title_size = 16,
                     plot_title_margin = 10,
                     ...) {
  
  ret <- ggplot2::theme_minimal(base_family = "Comic Sans MS",
                                base_size = base_size, ...) +
    scale_fill_manual(values=friendly_cols) +
    scale_color_manual(values=friendly_cols) +
    theme(panel.border = element_rect(fill = "NA",colour = "grey80")
    )
  ret$strip.text <- ggplot2::element_text(
    hjust = 0, size = strip_text_size,
    # face = 'bold',
    margin = ggplot2::margin(b = strip_text_margin),
    family = "Comic Sans MS"
  )
  ret$plot.subtitle <- ggplot2::element_text(
    hjust = 0, size = subtitle_size,
    margin = ggplot2::margin(b = subtitle_margin),
    family = "Comic Sans MS"
  )
  ret$plot.title <- ggplot2::element_text(
    hjust = 0, size = plot_title_size,
    margin = ggplot2::margin(b = plot_title_margin),
    family = "Comic Sans MS"
  ) 
  ret
  
}
  
  # Set theme
# my_theme <-
#   list(
#     scale_fill_manual(values=friendly_cols),
#     scale_color_manual(values=friendly_cols),
#     theme_bw() +
#       theme(
#         panel.border=element_blank(),
#         axis.line=element_line(),
#         panel.grid.major=element_line(size=0.2),
#         panel.grid.minor=element_line(size=0.1),
#         text=element_text(size=12),
#         legend.position="bottom",
#         aspect.ratio=1,
#         strip.background=element_blank(),
#         axis.title.x=element_text(margin=margin(t=10, r=10, b=10, l=10)),
#         axis.title.y=element_text(margin=margin(t=10, r=10, b=10, l=10))
#       )
#   )



# function to utilize the space from empty facet panels.
# https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2

library(gtable)
library(cowplot)

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}
