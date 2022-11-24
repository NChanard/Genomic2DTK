#' Create a viridis palette
#'
#' viridis
#' @keywords internal
#' @description Create a viridis palette.
#' @param paletteLength.num <numeric>: color number.
#' @param space.chr <numeric>: a character string; interpolation in RGB or CIE Lab color spaces. See ?grDevices::colorRamp for more details. (Default "rgb")
#' @param interpolate.chr <numeric>: use spline or linear interpolation. See ?grDevices::colorRamp for more details. (Default "linear")
#' @param bias.num <numeric>: a positive number.  Higher values give more widely spaced colors at the high end. See ?grDevices::colorRamp for more details. (Default 1)
#' @return  A vector of color.
#' @examples
#' viridis(9)
viridis <- function(paletteLength.num=NULL, space.chr='rgb', interpolate.chr='linear', bias.num=1){
  grDevices::colorRampPalette(
    colors=c("#440154", "#440256", "#450457", "#450559", "#46075A", "#46085C", "#460A5D",
      "#460B5E", "#470D60", "#470E61", "#471063", "#471164", "#471365", "#481467",
      "#481668", "#481769", "#48186A", "#481A6C", "#481B6D", "#481C6E", "#481D6F",
      "#481F70", "#482071", "#482173", "#482374", "#482475", "#482576", "#482677",
      "#482878", "#482979", "#472A7A", "#472C7A", "#472D7B", "#472E7C", "#472F7D",
      "#46307E", "#46327E", "#46337F", "#463480", "#453581", "#453781", "#453882",
      "#443983", "#443A83", "#443B84", "#433D84", "#433E85", "#423F85", "#424086",
      "#424186", "#414287", "#414487", "#404588", "#404688", "#3F4788", "#3F4889",
      "#3E4989", "#3E4A89", "#3E4C8A", "#3D4D8A", "#3D4E8A", "#3C4F8A", "#3C508B",
      "#3B518B", "#3B528B", "#3A538B", "#3A548C", "#39558C", "#39568C", "#38588C",
      "#38598C", "#375A8C", "#375B8D", "#365C8D", "#365D8D", "#355E8D", "#355F8D",
      "#34608D", "#34618D", "#33628D", "#33638D", "#32648E", "#32658E", "#31668E",
      "#31678E", "#31688E", "#30698E", "#306A8E", "#2F6B8E", "#2F6C8E", "#2E6D8E",
      "#2E6E8E", "#2E6F8E", "#2D708E", "#2D718E", "#2C718E", "#2C728E", "#2C738E",
      "#2B748E", "#2B758E", "#2A768E", "#2A778E", "#2A788E", "#29798E", "#297A8E",
      "#297B8E", "#287C8E", "#287D8E", "#277E8E", "#277F8E", "#27808E", "#26818E",
      "#26828E", "#26828E", "#25838E", "#25848E", "#25858E", "#24868E", "#24878E",
      "#23888E", "#23898E", "#238A8D", "#228B8D", "#228C8D", "#228D8D", "#218E8D",
      "#218F8D", "#21908D", "#21918C", "#20928C", "#20928C", "#20938C", "#1F948C",
      "#1F958B", "#1F968B", "#1F978B", "#1F988B", "#1F998A", "#1F9A8A", "#1E9B8A",
      "#1E9C89", "#1E9D89", "#1F9E89", "#1F9F88", "#1FA088", "#1FA188", "#1FA187",
      "#1FA287", "#20A386", "#20A486", "#21A585", "#21A685", "#22A785", "#22A884",
      "#23A983", "#24AA83", "#25AB82", "#25AC82", "#26AD81", "#27AD81", "#28AE80",
      "#29AF7F", "#2AB07F", "#2CB17E", "#2DB27D", "#2EB37C", "#2FB47C", "#31B57B",
      "#32B67A", "#34B679", "#35B779", "#37B878", "#38B977", "#3ABA76", "#3BBB75",
      "#3DBC74", "#3FBC73", "#40BD72", "#42BE71", "#44BF70", "#46C06F", "#48C16E",
      "#4AC16D", "#4CC26C", "#4EC36B", "#50C46A", "#52C569", "#54C568", "#56C667",
      "#58C765", "#5AC864", "#5CC863", "#5EC962", "#60CA60", "#63CB5F", "#65CB5E",
      "#67CC5C", "#69CD5B", "#6CCD5A", "#6ECE58", "#70CF57", "#73D056", "#75D054",
      "#77D153", "#7AD151", "#7CD250", "#7FD34E", "#81D34D", "#84D44B", "#86D549",
      "#89D548", "#8BD646", "#8ED645", "#90D743", "#93D741", "#95D840", "#98D83E",
      "#9BD93C", "#9DD93B", "#A0DA39", "#A2DA37", "#A5DB36", "#A8DB34", "#AADC32",
      "#ADDC30", "#B0DD2F", "#B2DD2D", "#B5DE2B", "#B8DE29", "#BADE28", "#BDDF26",
      "#C0DF25", "#C2DF23", "#C5E021", "#C8E020", "#CAE11F", "#CDE11D", "#D0E11C",
      "#D2E21B", "#D5E21A", "#D8E219", "#DAE319", "#DDE318", "#DFE318", "#E2E418",
      "#E5E419", "#E7E419", "#EAE51A", "#ECE51B", "#EFE51C", "#F1E51D", "#F4E61E",
      "#F6E620", "#F8E621", "#FBE723", "#FDE725"),
    space=space.chr,
    interpolate=interpolate.chr,
    bias=bias.num
    )(paletteLength.num)
}
