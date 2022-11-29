#' Create an Hue palette.
#'
#' Hue
#' @description Create an Hue palette.
#' @param paletteLength.num <numeric>: color number.
#' @param rotation.num <numeric>: if positive, rotates clockwise in the color space, reversing if the number is negative. If is NULL compute rotation according to hueRange.num parameter. (Default NULL)
#' @param hueRange.num <numeric>: Degree range in color space between 0 and 360. (Default c(0,360))
#' @param saturation.num <numeric>: Saturation value between 0 and 1. (Default 0.65)
#' @param lightness.num <numeric>: Lightness value between 0 and 1. (Default 0.65)
#' @param alpha.num <numeric>: Opacity value between 0 and 1. (Default 1)
#' @param alpha.bln <logical>: whether the alpha layer should be returned. (Default FALSE)
#' @return  A vector of color.
#' @examples
#' Hue(paletteLength.num=9)

Hue <- function(paletteLength.num=9, rotation.num=NULL, hueRange.num=c(0,360), saturation.num=0.65, lightness.num=0.65, alpha.num=1, alpha.bln=FALSE){
    if(paletteLength.num>2){
        if(is.null(rotation.num)){
            if(abs(diff(hueRange.num))>=180){
                rotation.num <- sign(diff(hueRange.num))
            }else{
                rotation.num <- -sign(diff(hueRange.num))
            }
        }
        if(sign(diff(hueRange.num))!=sign(rotation.num)){
            hueRange.num[which.min(hueRange.num)] <- hueRange.num[which.min(hueRange.num)] + 360
        }
        distGap.num = seq(hueRange.num[1],hueRange.num[2],length.out=paletteLength.num) |> diff() |> mean() |> abs()
        gap.num = min(abs(diff(hueRange.num)%%360),abs(360-abs(diff(hueRange.num)%%360)))
        if(paletteLength.num>1 && gap.num<distGap.num){
            adjust.num = ((distGap.num-gap.num)/(paletteLength.num+1) * paletteLength.num)
            hueRange.num[which.max(hueRange.num)] <- hueRange.num[which.max(hueRange.num)] - adjust.num
        }
        hue.lst <- seq(hueRange.num[1],hueRange.num[2],length.out=paletteLength.num)%%360
    }else if(paletteLength.num==2){
        hue.lst <- Hue(paletteLength.num=5, rotation.num=rotation.num, hueRange.num=hueRange.num, saturation.num=saturation.num, lightness.num=lightness.num, alpha.num=alpha.num, alpha.bln=alpha.bln)[c(2,3)]
    }else if(paletteLength.num==1){
        hue.lst <- Hue(paletteLength.num=3, rotation.num=rotation.num, hueRange.num=hueRange.num, saturation.num=saturation.num, lightness.num=lightness.num, alpha.num=alpha.num, alpha.bln=alpha.bln)[2]
    }
    if(is.numeric(hue.lst)){
        hue.lst <- lapply(hue.lst,function(hue.num){Hsl2Hex(c(hue=hue.num,saturation=saturation.num,light=lightness.num,alpha=alpha.num), alpha.bln=alpha.bln)}) |> unlist()
        return(hue.lst)
    }else{
        return(hue.lst)
    }
}
