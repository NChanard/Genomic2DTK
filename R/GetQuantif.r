#' GetQuantif
#'
#' Compute a function on each matrix that compose a matrices list.
#' @param matrices.lst <List[matrix]> : A matrices list.
#' @param area.fun <character or function> : A character or function thaht alow an extraction of an area on each matrix that compose the matrices lits (Default "center").
#' \itemize{
#' \item "mean_rm0" : apply a mean after replace 0 with NA
#' \item "median_rm0" : apply a median after replace 0 with NA
#' \item "sum_rm0" : apply a sum after replace 0 with NA
#' \item "median" : apply a median
#' \item "sum" : apply a sum
#' \item "mean" or other character :  apply a mean
#' }
#' @param operation.fun <character or function> : A character or function to applt on the  extracted area of each matrix that compose the matrices lits (Default "mean_rm0").
#' \itemize{
#' \item "C" or "CENTER" : Take the area of interaction between anchor and bait.
#' \item "UL" or "UPPER_LEFT" : Take the area of interaction in the uppper left square
#' \item "UR" or "UPPER_RIGHT" : Take the area of interaction in the uppper right square
#' \item "BL" or "BOTTOM_LEFT" : Take the area of interaction in the bottom left square
#' \item "BR" or "BOTTOM_RIGHT" : Take the area of interaction in the bottom right square
#' \item "U" or "UPPER" : Take the area of interaction in above the center area
#' \item "B" or "BOTTOM" : Take the area of interaction in below the center area
#' \item "L" or "LEFT" : Take the area of interaction in the left of the center area
#' \item "R" or "RIGHT" : Take the area of interaction in the right of the center area
#' \item "D" or "DONUT" : Take the area of interaction that surround the center area
#' }
#' @param name.chr <character> : The name of a column in GInteraction attributes of matrices.lst used as named in the output (Default NULL).
#' @return A GRange object.

GetQuantif = function(matrices.lst, area.fun="center", operation.fun="mean_rm0", name.chr=NULL){
    # Define operation function
        if(is.null(operation.fun)){
            operation.fun <- function(area){c(area[which(!is.na(area))])}
        }else if(!is.function(operation.fun)){
            operation.fun <- dplyr::case_when(
                operation.fun == "mean_rm0"                     ~ "function(x){x[which(x==0)]<-NA;mean(x,na.rm=TRUE)}",
                operation.fun == "median_rm0"                   ~ "function(x){x[which(x==0)]<-NA;stats::median(x,na.rm=TRUE)}",
                operation.fun == "sum_rm0"                      ~ "function(x){x[which(x==0)]<-NA;sum(x,na.rm=TRUE)}",
                operation.fun == "median"                       ~ "function(x){stats::median(x,na.rm=TRUE)}",
                operation.fun == "sum"                          ~ "function(x){sum(x,na.rm=TRUE)}",
                operation.fun == "mean"                         ~ "function(x){mean(x,na.rm=TRUE)}",
                TRUE                                            ~ "function(x){mean(x,na.rm=TRUE)}"
                ) %>% DataTK::WrapFunction(.)
        }
    # Define extraction function
        matriceDim.num = attributes(matrices.lst)$matriceDim
        if(!is.function(area.fun) && attributes(matrices.lst)$referencePoint == "pf"){
            # Compute rows and cols index
                # Center
                    center.num = c(floor((matriceDim.num+1)/2) + ifelse(matriceDim.num>=9,yes=1,no=0) , ceiling((matriceDim.num+1)/2) - ifelse(matriceDim.num>=9,yes=1,no=0))
                    centerStart.num = min(center.num)
                    centerEnd.num = max(center.num)
                    center.chr = ifelse(centerStart.num==centerEnd.num,yes=centerStart.num,no=paste0( centerStart.num, ":", centerEnd.num ))
                # Rest
                    first.chr = paste0("1:", centerStart.num - 2)
                    second.chr = paste0(centerEnd.num + 2, ":", matriceDim.num)
            # Compute rows and cols index for Donut
                # Thick
                    donutThick.num = ifelse(matriceDim.num>=9,yes=2,no=1)
                # Top
                    topDonutStart.num = (centerStart.num - 1 - donutThick.num)
                    topDonutEnd.num = centerStart.num - 2
                    topDonut.num = topDonutStart.num:topDonutEnd.num
                # Left
                    leftDonut.num = topDonut.num
                # Bot
                    botDonutEnd.num = centerEnd.num + 1 + donutThick.num
                    botDonutStart.num = centerEnd.num + 2
                    botDonut.num = botDonutStart.num:botDonutEnd.num
                # Right
                    rightDonut.num = botDonut.num
                # Width
                    widthDonut.num = (centerStart.num - 1) : (centerEnd.num + 1)
                # Height
                    heigthDonut.num = widthDonut.num
                # Coords
                    donutCoord.dtf = rbind(
                        expand.grid(topDonut.num,leftDonut.num),
                        expand.grid(topDonut.num,widthDonut.num),
                        expand.grid(topDonut.num,rightDonut.num),
                        expand.grid(heigthDonut.num,leftDonut.num),
                        expand.grid(heigthDonut.num,rightDonut.num),
                        expand.grid(botDonut.num,leftDonut.num),
                        expand.grid(botDonut.num,widthDonut.num),
                        expand.grid(botDonut.num,rightDonut.num)
                    )
                    donutRows.chr =  paste(donutCoord.dtf[,1],collapse=",")
                    donutCols.chr =  paste(donutCoord.dtf[,2],collapse=",")
                    donut.chr = paste0("cbind(  c(", donutRows.chr,")  ,  c(",donutCols.chr,")  )")
            # Wrap index in a function
                area.fun <- dplyr::case_when(
                    toupper(area.fun) %in% c("C","CENTER")          ~ list(center.chr, center.chr),
                    toupper(area.fun) %in% c("UL","UPPER_LEFT")     ~ list(first.chr , first.chr ),
                    toupper(area.fun) %in% c("UR","UPPER_RIGHT")    ~ list(first.chr , second.chr),
                    toupper(area.fun) %in% c("BL","BOTTOM_LEFT")    ~ list(second.chr, first.chr ),
                    toupper(area.fun) %in% c("BR","BOTTOM_RIGHT")   ~ list(second.chr, second.chr),
                    toupper(area.fun) %in% c("U","UPPER")           ~ list(first.chr , center.chr),
                    toupper(area.fun) %in% c("B","BOTTOM")          ~ list(second.chr, center.chr),
                    toupper(area.fun) %in% c("L","LEFT")            ~ list(center.chr, first.chr ),
                    toupper(area.fun) %in% c("R","RIGHT")           ~ list(center.chr, second.chr),
                    toupper(area.fun) %in% c("D","DONUT")           ~ list(donut.chr             ),
                    TRUE                                            ~ list(center.chr, center.chr)
                ) %>% paste(collapse=",") %>% paste0("function(matrice.mtx){ matrice.mtx[",.,"] }" ) %>%  DataTK::WrapFunction(.)
        }else if(!is.function(area.fun) && attributes(matrices.lst)$referencePoint == "rf"){
            shiftFactor = attributes(matrices.lst)$shiftFactor
            # Compute rows and cols index
                # Anchor
                    anchorStart.num = max(1,floor((matriceDim.num-2)*shiftFactor/(1+2*shiftFactor))+1 - ifelse(matriceDim.num>=9,yes=1,no=0))
                    anchorEnd.num = max(1,min(matriceDim.num,floor((matriceDim.num-2)*shiftFactor/(1+2*shiftFactor))+ 1 + ifelse(matriceDim.num>=9,yes=1,no=0)))
                    anchor.chr = ifelse(anchorStart.num==anchorEnd.num, yes=anchorStart.num, no = paste0(anchorStart.num,":",anchorEnd.num))
                # Bait
                    baitStart.num = max(1,ceiling((matriceDim.num-2)*(1+shiftFactor)/(1+2*shiftFactor))+ 2 - ifelse(matriceDim.num>=9,yes=1,no=0))
                    baitEnd.num = max(1,min(matriceDim.num,ceiling((matriceDim.num-2)*(1+shiftFactor)/(1+2*shiftFactor))+ 2 + ifelse(matriceDim.num>=9,yes=1,no=0)))
                    bait.chr = ifelse(baitStart.num==baitEnd.num, yes=baitStart.num, no = paste0(baitStart.num,":",baitEnd.num))
                # Rest
                    first.chr = paste0("1:", anchorStart.num - 2)
                    second.chr = paste0(baitEnd.num + 2, ":", matriceDim.num)
                    ULwidth.chr = paste0("1:", baitStart.num - 2 )
                    inner.chr = paste0(anchorEnd.num  +  2,":", baitStart.num  - 2)
                    BRheight.chr = paste0(anchorEnd.num +  2 , ":", matriceDim.num)
            # Computability
                U.lgk <- anchorStart.num >= 3
                R.lgk <- matriceDim.num >= baitEnd.num+2
                B.lgk <- (((matriceDim.num+1)/2)-anchorEnd.num) >= 1
                L.lgk <- (baitStart.num-((matriceDim.num+1)/2)) >= 1
                UL.lgk <- U.lgk && baitStart.num >= 3
                UR.lgk <- U.lgk && R.lgk
                BR.lgk <- matriceDim.num >= anchorEnd.num+2 && R.lgk
                BL.lgk <- (((matriceDim.num+1)/2)-anchorEnd.num) >= 2 && (baitStart.num-((matriceDim.num+1)/2)) >= 2
                D.lgk <- sum(U.lgk, R.lgk, B.lgk, L.lgk, UL.lgk, UR.lgk, BR.lgk, BL.lgk) >=1
            # Compute rows and cols index for Donut
                if(D.lgk){
                    # Thick
                        donutThick.num = ifelse(matriceDim.num>=9,yes=2,no=1)
                    # Top
                        topDonutStart.num = (anchorStart.num - 1 - donutThick.num)
                        topDonutEnd.num = anchorStart.num - 2
                        topDonut.num = topDonutStart.num:topDonutEnd.num
                    # Left
                        leftDonutStart.num = (baitStart.num - 1 - donutThick.num)
                        leftDonutEnd.num = baitStart.num - 2
                        leftDonut.num = leftDonutStart.num:leftDonutEnd.num
                    # Bot
                        botDonutEnd.num = anchorEnd.num + 1 + donutThick.num
                        botDonutStart.num = anchorEnd.num + 2
                        botDonut.num = botDonutStart.num:botDonutEnd.num
                    # Right
                        rightDonutEnd.num = baitEnd.num + 1 + donutThick.num
                        rightDonutStart.num = baitEnd.num + 2
                        rightDonut.num = rightDonutStart.num:rightDonutEnd.num
                    # Widht
                        widthDonut.num = (baitStart.num - 1) : (baitEnd.num + 1)
                    # Height
                        heigthDonut.num = (anchorStart.num - 1) : (anchorEnd.num + 1)
                    # Coord
                        donutCoord.dtf = rbind(
                                expand.grid(topDonut.num,leftDonut.num),
                                expand.grid(topDonut.num,widthDonut.num),
                                expand.grid(topDonut.num,rightDonut.num),
                                expand.grid(heigthDonut.num,leftDonut.num),
                                expand.grid(heigthDonut.num,rightDonut.num),
                                expand.grid(botDonut.num,leftDonut.num),
                                expand.grid(botDonut.num,widthDonut.num),
                                expand.grid(botDonut.num,rightDonut.num)
                            ) %>% dplyr::filter(Var1>=1 & Var2<=matriceDim.num & Var1<=Var2)
                        donutRows.chr =  paste(donutCoord.dtf[,1],collapse=",")
                        donutCols.chr =  paste(donutCoord.dtf[,2],collapse=",")
                        donut.chr = paste0("cbind(  c(", donutRows.chr,")  ,  c(",donutCols.chr,")  )")
                }
            # Wrap index in a function
                area.fun <- dplyr::case_when(
                    toupper(area.fun) %in% c("C","CENTER")                  ~ list(anchor.chr   , bait.chr     ),
                    toupper(area.fun) %in% c("UL","UPPER_LEFT")   && UL.lgk ~ list(first.chr , ULwidth.chr  ),
                    toupper(area.fun) %in% c("UR","UPPER_RIGHT")  && UR.lgk ~ list(first.chr , second.chr),
                    toupper(area.fun) %in% c("BL","BOTTOM_LEFT")  && BL.lgk ~ list(inner.chr, inner.chr        ),
                    toupper(area.fun) %in% c("BR","BOTTOM_RIGHT") && BR.lgk ~ list(BRheight.chr, second.chr),
                    toupper(area.fun) %in% c("U","UPPER")         && U.lgk  ~ list(first.chr    , bait.chr   ),
                    toupper(area.fun) %in% c("B","BOTTOM")        && B.lgk  ~ list(second.chr   , bait.chr   ),
                    toupper(area.fun) %in% c("L","LEFT")          && L.lgk  ~ list(anchor.chr   , first.chr    ),
                    toupper(area.fun) %in% c("R","RIGHT")         && R.lgk  ~ list(anchor.chr   , second.chr   ),
                    toupper(area.fun) %in% c("D","DONUT")         && D.lgk  ~ list(donut.chr                   ),
                    TRUE                                                    ~ list(anchor.chr   , bait.chr     )
                ) %>% paste(collapse=",") %>% paste0("function(matrice.mtx){ matrice.mtx[",.,"] }" ) %>%  DataTK::WrapFunction(.)   
        }
    # Compute quantif
        quantif.num = lapply(matrices.lst, function(mtx){ 
            mtxQuantif.num <- mtx %>%
                area.fun(.) %>%
                operation.fun(.) %>%
                magrittr::set_names(., NULL)
            rownames(mtxQuantif.num) <- NULL
            colnames(mtxQuantif.num) <- NULL
            return(mtxQuantif.num)
        })
    # Get Names
        if(!is.null(name.chr)){
            names.chr_lst <- attributes(matrices.lst)$interactions %>%
                S4Vectors::mcols(.) %>%
                data.frame %>%
                dplyr::arrange(factor(submatrix.name, levels=names(quantif.num))) %>% dplyr::pull(name.chr)
        } else {
            names.chr_lst <- names(quantif.num)
        }
    # Repeted Index if names.chr_lst is a nested List
        repeted.ndx <- names.chr_lst %>% lapply(.,length) %>% magrittr::set_names(seq_along(.)) %>% rep(names(.),.) %>% as.numeric
    # Add attributes
        quantif.num[repeted.ndx] %>%
            magrittr::set_names(unlist(names.chr_lst )) %>% 
            DevTK::AddAttr(c(
                attributes(matrices.lst),
                interactions = attributes(matrices.lst)$interactions[repeted.ndx],
                operation = operation.fun,
                area = area.fun,
                duplicated = which(duplicated(repeted.ndx))
            )) %>% return(.y)
}