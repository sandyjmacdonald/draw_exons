import random
from PIL import Image, ImageDraw
from draw_exons import Exon, CDS, CDS_group

# Example showing how to draw multiple CDSs in a group,
# displaying each in a different colour.

# Randomly pick a number of CDSs and total exons for the gene.
num_cdss = random.randint(7, 12)
num_exons = random.randint(7, 15)

# Empty list to store the gene in.
gene = []

# Empty list to store exon coords for each CDS.
cds_coords = []

# Keep track of the current position when laying
# out exons.
curr_pos = 0

# Colorbrewer palette of 12 qualitative colours.
pal = [
    (141, 211, 199),
    (255, 255, 179),
    (190, 186, 218),
    (251, 128, 114),
    (128, 177, 211),
    (253, 180, 98),
    (179, 222, 105),
    (252, 205, 229),
    (217, 217, 217),
    (188, 128, 189),
    (204, 235, 197),
    (255, 237, 111),
]

# Build exons and add to gene list.
for j in range(num_exons):
    exon_length = random.randint(100, 1000)
    gap = random.randint(100, 1000)
    start = curr_pos + gap
    end = start + exon_length
    curr_pos = end
    gene.append((start, end))

# For each CDS, randomly pick some exons and add
# to the cds_coords list.
for i in range(num_cdss):
    chosen_exons = sorted(
        random.sample(gene, random.randint(2, num_exons)), key=lambda x: x[0]
    )
    cds_coords.append(chosen_exons)

# Randomly pick some colours from the palette.
fill_colours = random.sample(pal, len(cds_coords))

# Just use black for the outline colours.
outline_colours = [(0, 0, 0) for i in range(len(cds_coords))]

# Scale up size, and then resample down later to improve
# jaggy lines.
scale_factor = 4
cds_height = 20 * scale_factor
width = 800 * scale_factor

# Create a CDS_group with the cds_coords and colours.
CDS_group = CDS_group(
    cds_coords=cds_coords,
    cds_height=cds_height,
    cds_spacing=int(cds_height * 1.5),
    linewidth=5,
    width=width,
    fill_colours=fill_colours,
    outline_colours=outline_colours,
)

# Draw the CDSs
img = CDS_group.draw()

# Resample and size down to improve jaggy lines.
w, h = img.size
img_1x = img.resize(
    (w // scale_factor, h // scale_factor), resample=Image.Resampling.BICUBIC
)

# Show the CDSs.
img_1x.show()
