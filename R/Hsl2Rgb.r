#' Convert a color in HSL to RGB.
#'
#' Hsl2Rgb
#' @keywords internal
#' @description Convert a color in HSl (Hue,Saturation,Light) format to RGB format.
#' @param hsl.col <charcater>: a vector of the color's HSL code.
#' @param alpha.bln <logical>: whether the alpha layer should be returned. (Default FALSE)
#' @return An integer vector of the color's RGB code.
#' @examples
#' Hsl2Rgb(c(43.8,0.873,0.492,0.498),alpha.bln=TRUE)
Hsl2Rgb <- function(hsl.col=NULL, alpha.bln=FALSE){
    if(3>length(hsl.col) | length(hsl.col)>4){
        cat("Error, need 3 or 4 values beetween 0 and 255, first value for hue, second for saturation, third for light and last for alpha")
    }else{
        if(IsHsl(hsl.col)){
            if(length(hsl.col) == 3){
                alpha.num <- 255
            }else{
                alpha.num <- hsl.col[4]*255
                hsl.col <- hsl.col[1:3]
            }
            C <- hsl.col[2] * (1-abs(2*hsl.col[3]-1))
            X <- C * (1-abs((hsl.col[1]/60)%%2-1))
            m <- hsl.col[3] - C/2
            if (0<=hsl.col[1] & hsl.col[1]<60){
                rgb.col <- c( red=(C+m)*255, green=(X+m)*255, blue=(0+m)*255 )
            }else if(60<=hsl.col[1] & hsl.col[1]<120){
                rgb.col <- c( red=(X+m)*255, green=(C+m)*255, blue=(0+m)*255 )
            }else if(120<=hsl.col[1] & hsl.col[1]<180){
                rgb.col <- c( red=(0+m)*255, green=(C+m)*255, blue=(X+m)*255 )
            }else if(180<=hsl.col[1] & hsl.col[1]<240){
                rgb.col <- c( red=(0+m)*255, green=(X+m)*255, blue=(C+m)*255 )
            }else if(240<=hsl.col[1] & hsl.col[1]<300){
                rgb.col <- c( red=(X+m)*255, green=(0+m)*255, blue=(C+m)*255 )
            }else if(300<=hsl.col[1] & hsl.col[1]<360){
                rgb.col <- c( red=(C+m)*255, green=(0+m)*255, blue=(X+m)*255 )
            }
            if(alpha.bln){ rgb.col <- c(rgb.col, alpha=alpha.num)}
        }else{
            cat("Error, need 3 or 4 values beetween 0 and 255, first value for hue, second for saturation, third for light and last for alpha")
        }
    }
    return(round(rgb.col))
}
        