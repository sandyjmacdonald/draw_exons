from PIL import Image, ImageDraw
from draw_exons import Exon, CDS, CDS_group, add_margin

# Example showing how to draw several CDSs and add descriptive
# labels alongside them.

# Exon start/end positions in gene.
gene = [(80, 420), (600, 770), (900, 1220), (1500, 1730), (1800, 1900)]

# Palette of 4 qualitative colours.
fill_colours = [(141, 211, 199), (255, 255, 179), (190, 186, 218), (251, 128, 114)]

# Dictionary of the transcript labels and exons.
cdss = {
    "Transcript A": [1, 2, 3, 4, 5],
    "Transcript B": [1, 3, 4],
    "Transcript C": [2, 3, 4, 5],
    "Transcript D": [1, 4, 5],
}

# Empty list to add the CDS coordinates to.
cds_coords = []

# Empty list to add the transcript labels to.
cds_ids = []

# Go through the transcripts and add the labels and exon
# coordinates to their respective lists.
for label, coords in cdss.items():
    cds_ids.append(label)
    exon_coords = [gene[i - 1] for i in coords]
    cds_coords.append(exon_coords)

# Just use black for the outline colours.
outline_colours = [(0, 0, 0) for i in range(len(cds_coords))]

# Scale up size, and then resample down later to improve
# jaggy lines.
scale_factor = 4
cds_height = 20 * scale_factor
width = 800 * scale_factor

# Size (in points) of labels for transcripts.
label_size = 22
label_size *= scale_factor

# Add some margin for the labels and trim down again after.
margin = 150 * scale_factor
trim = -100 * scale_factor

# Create a CDS_group with the cds_coords, labels, and colours.
CDS_group = CDS_group(
    cds_coords=cds_coords,
    cds_height=cds_height,
    cds_ids=cds_ids,
    cds_spacing=int(cds_height * 1.5),
    linewidth=5,
    width=width,
    fill_colours=fill_colours,
    outline_colours=outline_colours,
    margin=margin,
)

# Draw the CDSs and trim extra margin.
img = CDS_group.draw(label=True, label_size=label_size)
img = add_margin(img, top=trim, right=trim, bottom=trim, left=trim)

# Resample and size down to improve jaggy lines.
w, h = img.size
img_1x = img.resize(
    (w // scale_factor, h // scale_factor), resample=Image.Resampling.BICUBIC
)

# Show the CDSs.
img_1x.show()
