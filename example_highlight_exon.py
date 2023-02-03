import random
from PIL import Image, ImageDraw
from draw_exons import Exon, CDS, CDS_group

# A simple example illustrating how to draw a single CDS
# with one exon highlighted in a different colour to the rest.

# Scale image up and then scale down to fix jaggy lines.
scale_factor = 4
w, h = (x * scale_factor for x in (800, 200))

# Create a new image and canvas to draw the CDS on.
img = Image.new("RGB", (w, h), "white")
canv = ImageDraw.Draw(img, mode="RGBA")

# The start/end coordinates of each exon.
exon_coords = [(80, 420), (600, 770), (900, 1220), (1500, 1730), (1800, 1900)]

# Main fill colour for the exons. They'll all be this
# colour initially.
fill = (141, 211, 199)

linewidth = 2

# Create the CDS
cds = CDS(
    canv=canv,
    exon_coords=exon_coords,
    x1=100 * scale_factor,
    x2=700 * scale_factor,
    y=75 * scale_factor,
    height=50 * scale_factor,
    fill=fill,
    linewidth=linewidth * scale_factor,
)

# Draw the CDS.
cds.draw()

# Change the fill colour of exon 3 (index 2), and
# redraw the CDS.
cds.exons[2].fill = (190, 186, 218)
cds.redraw()

# Scale the image back down, and show it.
img_1x = img.resize(
    (w // scale_factor, h // scale_factor), resample=Image.Resampling.BICUBIC
)
img_1x.show()
