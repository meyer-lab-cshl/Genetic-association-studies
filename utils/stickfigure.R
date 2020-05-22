draw.stick <- function(x,y,scale=1,arms="down", 
                       gender="male",lwd=3, clcol="white",
                       face="happy", linecol=gray(.3),
                       hat=NA) {
  # clcol: color of clothes - any color
  # scale: fize of figure
  # x,y left bottom alignment of figure
  # linecol: color of lines - any color
  # lwd: line weight
  # Arms: "down", "nuetral", "up", "hip", "wave"
  # Gender: "male", "female"
  # Face: "happy", "sad", "annoyed", "surprised"
  # Hat: plot hat T,F


  # Set the figure scale, default it 1
  s <- scale/100

  # If is undefined then give the man a hat
  if (is.na(hat)) hat<-(gender=="male")

  require("plotrix")

  # Draw Head
  draw.ellipse(x+50*s,y+75*s,10*s,15*s,lwd=lwd, border=linecol)

  if (face=="happy") {
    # Draw eyes
    draw.ellipse(x+46*s,y+77*s,2.5*s,2*s,lwd=lwd, border=linecol)
    draw.ellipse(x+54*s,y+77*s,2.5*s,2*s,lwd=lwd, border=linecol)

    # Draw mouth
    draw.ellipse(x+50*s,y+72*s,6*s,8*s, segment = c(-160,-20),
                  lwd=lwd, border=linecol)
  }
  if (face=="sad")   {
    # Draw eyes
    draw.ellipse(x+46*s,y+75*s,2*s,2*s,lwd=lwd, border=linecol)
    draw.ellipse(x+54*s,y+75*s,2*s,2*s,lwd=lwd, border=linecol)

    # Draw mouth
    draw.ellipse(x+50*s,y+60*s,6*s,8*s, segment = c(140,40),
                 lwd=lwd, border=linecol)
  }

  if (face=="surprised")   {
    # Draw eyes
    draw.ellipse(x+46*s,y+78*s,3*s,2*s,lwd=lwd, border=linecol)
    draw.ellipse(x+54*s,y+78*s,3*s,2*s,lwd=lwd, border=linecol)

    # Draw irises
  #  draw.ellipse(x+46*s,y+78*s,1*s,1*s,lwd=lwd, border=linecol)
  #  draw.ellipse(x+54*s,y+78*s,1*s,1*s,lwd=lwd, border=linecol)

    # Draw mouth
    draw.ellipse(x+50*s,y+65*s,3*s,4*s,
                 lwd=lwd, border=linecol)
  }


  if (face=="annoyed")  {
    # Draw mouth
    lines(c(x+46*s,x+55*s), c(y+66*s,y+68*s),lwd=lwd, col=linecol)

    # Draw eyes
    draw.ellipse(x+46*s,y+76*s,2*s,2*s,lwd=lwd, border=linecol)
    draw.ellipse(x+54*s,y+76*s,2*s,1*s,lwd=lwd, border=linecol)
  }




  # Draw torso
  lines(c(x+50*s,x+50*s), c(y+35*s,y+60*s),lwd=lwd, col=linecol)

  # Draw arms
  if (arms=="down") {
    lines(c(x+50*s,x+36*s), c(y+55*s,y+30*s),lwd=lwd, col=linecol) # Left
    lines(c(x+50*s,x+64*s), c(y+55*s,y+30*s),lwd=lwd, col=linecol) # Right
  }
  if (arms=="nuetral") {
    lines(c(x+50*s,x+30*s), c(y+50*s,y+55*s),lwd=lwd, col=linecol) # Left
    lines(c(x+50*s,x+70*s), c(y+50*s,y+55*s),lwd=lwd, col=linecol) # Right
  }
  if (arms=="up") {
    lines(c(x+50*s,x+32*s), c(y+50*s,y+65*s),lwd=lwd, col=linecol) # Left
    lines(c(x+50*s,x+68*s), c(y+50*s,y+65*s),lwd=lwd, col=linecol) # Right
  }
  if (arms=="hip") {
    lines(c(x+50*s,x+37*s,x+48*s), c(y+56*s,y+47*s,y+40*s),lwd=lwd, col=linecol) # Left
    lines(c(x+50*s,x+63*s,x+51*s), c(y+56*s,y+49*s,y+62*s),lwd=lwd, col=linecol) # Right
  }
  if (arms=="wave") {
    lines(c(x+50*s,x+38*s,x+33*s), c(y+56*s,y+60*s,y+78*s),lwd=lwd, col=linecol) # Left
    lines(c(x+50*s,x+63*s,x+52*s), c(y+56*s,y+47*s,y+40*s),lwd=lwd, col=linecol) # Right
  }

  # Draw male legs
  if (gender=="male") {
    lines(c(x+50*s,x+40*s), c(y+35*s,y+5*s),lwd=lwd, col=linecol)
    lines(c(x+50*s,x+60*s), c(y+35*s,y+5*s),lwd=lwd, col=linecol)

  }

  # Draw female legs and dress
  if (gender!="male") {
    # Draw legs
    lines(c(x+45*s,x+45*s), c(y+17*s,y+5*s),lwd=lwd, col=linecol)
    lines(c(x+55*s,x+55*s), c(y+17*s,y+5*s),lwd=lwd, col=linecol)

    # Draw dress
    polygon(c(x+s*50,x+s*35,x+s*65), 
          c(y+s*40,y+s*17,y+s*17), 
            col=clcol, border=linecol,lwd=lwd)
  }

  # Draw hat
  if (hat==T) polygon(
    c(x+35*s,x+65*s,x+65*s,x+59*s,x+58*s,x+42*s,x+41*s,x+35*s),
    c(y+84*s,y+84*s,y+86*s,y+86*s,y+91*s,y+91*s,y+86*s,y+86*s),
    col=clcol, border=linecol,lwd=lwd)

}

# Map of stick figures
par(mar = rep(0, 4))

plot(c(.25,1.25), c(0,3), type="n", xaxt='n', yaxt='n', ann=FALSE)

  draw.stick(0,0, arms="hip")
  draw.stick(.5,0, gender="female", arms="up")

  draw.stick(0,1, arms="nuetral", lwd=2, linecol=gray(.5), 
             clcol="red", face="sad")

  draw.stick(.5,1, gender="female", arms="down", clcol="purple", 
             lwd=2, linecol=gray(.5), face="sad",hat=T)

  draw.stick(0,2, arms="wave", linecol=gray(.7),
             clcol="blue", face="surprised",hat=F)

  draw.stick(.5,2, gender="female", arms="hip", clcol="light blue", 
             linecol=gray(.7),face="annoyed")

