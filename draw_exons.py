import random
from PIL import Image, ImageDraw


class Exon:
    """
    Class representing a single exon.

    Attributes
    ----------
    canv: a PIL ImageDraw.Draw() canvas object
        the canvas on which the exon is drawn
    x1=: int
        the left-most x-coordinate of the exon
    x2=: int
        the right-most x-coordinate of the exon
    y=: int
        the y-coordinate of the centre line of the exon
    height=: int
        the height in pixels of the exon
    fill=: 3-tuple of 0-255
    """

    def __init__(
        self,
        canv,
        x1=0,
        x2=100,
        y=0,
        height=10,
        fill=(255, 255, 255),
        outline=(0, 0, 0),
        linewidth=1,
    ):
        self.canv = canv
        self.x1 = x1
        self.x2 = x2
        self.y = y
        self.height = height
        self.fill = fill
        self.outline = outline
        self.linewidth = linewidth

    def draw(self):
        """
        Draws the exon on a PIL canvas
        """
        self.canv.rectangle(
            (self.x1, self.y, self.x2, self.y + self.height),
            fill=self.fill,
            outline=self.outline,
            width=self.linewidth,
        )


class CDS:
    """
    Class representing a group of exons in a coding sequence.

    Attributes
    ----------
    canv: a PIL ImageDraw.Draw() canvas object
        the canvas on which the CDS is drawn
    cds_id=: str
        an ID for the CDS
    exon_coords=: list of 2-tuples
        a list of the start and end coordinates
        for all of the exons in CDS, as tuples
    x1=: int
        the left-most x-coordinate of the CDS
    x2=: int
        the right-most x-coordinate of the CDS
    y=: int
        the y-coordinate of the centre line of the CDS
    height=: int
        the height in pixels of the CDS
    fill=: 3-tuple of 0-255
        the fill colour for the CDS
    outline=: 3-tuple of 0-255
        the outline colour for the CDS
    linewidth=: int
        the width of the outline of the CDS
    """

    def __init__(
        self,
        canv,
        cds_id=None,
        exon_coords=[],
        x1=0,
        x2=100,
        y=0,
        height=10,
        fill=(255, 255, 255),
        outline=(0, 0, 0),
        linewidth=1,
    ):
        self.canv = canv
        self.cds_id = cds_id
        self.exon_coords = exon_coords
        self.x1 = x1
        self.x2 = x2
        self.y = y
        self.height = height
        self.width = self.x2 - self.x1
        self.fill = fill
        self.outline = outline
        self.linewidth = linewidth
        self.exons = []

    def make_exons(self):
        """
        Returns a list of Exon instances for the CDS.
        """
        all_coords = [coord for coords in self.exon_coords for coord in coords]
        min_x = min(all_coords)
        max_x = max(all_coords)
        range_x = max_x - min_x
        for coords in self.exon_coords:
            start, end = coords
            x1 = self.x1 + int(((start - min_x) / range_x) * self.width)
            x2 = self.x1 + int(((end - min_x) / range_x) * self.width)
            exon = Exon(
                self.canv,
                x1=x1,
                x2=x2,
                y=self.y,
                height=self.height,
                fill=self.fill,
                outline=self.outline,
                linewidth=self.linewidth,
            )
            self.exons.append(exon)
        return self.exons

    def draw_exons(self):
        """
        Calls the draw method of the Exon class instances
        and returns the resulting canvas object for further
        processing.
        """
        for exon in self.exons:
            exon.draw()
        return self.canv

    def draw_lines(self):
        """
        Draws lines (introns) between the exons in the CDS,
        and returns the modified PIL canvas.
        """
        for i in range(1, len(self.exons)):
            last_exon = self.exons[i - 1]
            this_exon = self.exons[i]
            x1 = last_exon.x2
            x2 = this_exon.x1
            midpoint = int(x1 + ((x2 - x1) / 2))
            y = int(self.y + (self.height / 2))
            y_midpoint = y - (self.height / 4)
            linewidth = self.linewidth
            self.canv.line(
                (x1, y, midpoint, y_midpoint), fill=self.outline, width=linewidth
            )
            self.canv.line(
                (midpoint, y_midpoint, x2, y), fill=self.outline, width=linewidth
            )
        return self.canv

    def draw_label(self, x, y, label_size=16, colour=(0, 0, 0), centre=True):
        """
        Add a label for the CDS, centred on the y-coordinate.
        Returns the coordinates of the bounding box for the
        text label.

            Arguments:
                x (int): pixel coordinate for the left hand edge of the label
                y (int): pixel coordinate for the vertical alignment of the label

            Parameters:
                label_size (int): the font size in points for the label (default 16)
                colour (3-tuple of ints): 0-255 values for RGB colour of the text (default (255, 255, 255))
                centre (Boolean): whether to centre the label vertically

            Returns:
                bounding_box (4-tuple of tuples): the bounding box of the text

        """
        from PIL import ImageFont

        font = ImageFont.truetype("Roboto-Regular.ttf", label_size)
        left, top, right, bottom = font.getbbox(self.cds_id)
        if centre:
            label_height = bottom - top
            cds_height = self.height - (0.5 * self.linewidth)
            offset = int((cds_height - label_height) / 2)
            y = int(y + offset - (0.5 * self.linewidth))
        self.canv.text((x, y), self.cds_id, colour, font=font)

        bounding_box = (left, top, right, bottom)

        return bounding_box

    def add_exon(self, exon):
        """
        Add an exon to the CDS.

            Arguments:
                exon (Exon() instance): an Exon() class instance
        """
        self.exons.append(exon)

    def draw(self):
        """
        A do-all function that creates the exon instances,
        draws them, then draw the connecting lines. Returns
        the modified canvas object.
        """
        self.make_exons()
        self.draw_exons()
        self.draw_lines()

        return self.canv

    def redraw(self):
        """
        A do-all function that creates the exon instances,
        draws them, then draw the connecting lines. Returns
        the modified canvas object.
        """
        this_canv = self.canv
        w, h = this_canv._image.size
        new_img = Image.new("RGB", (w, h), "white")
        new_canv = ImageDraw.Draw(new_img, mode="RGBA")
        self.canv = new_canv
        self.make_exons()
        self.draw_exons()
        self.draw_lines()

        return self.canv


class CDS_group:
    """
    Class representing a group of CDS() instances.

    Attributes
    ----------
    cds_coords=: list
        list of lists of exon start/end coordinates for each CDS
    cds_ids=: list
        list of strings with IDs for each CDS
    cdss=: list
        list of CDS() instances for the CDSs to be drawn
    cds_height=: int
        the height in pixels of each CDS (default 50)
    cds_spacing=: int
        the padding between each CDS vertically in pixels (default 50)
    width=: int
        the width of the whole CDS group in pixels (default 800)
    margin=: int
        the size in pixels of the margin in pixels (default 50)
    fill_colours=: list
        3-tuples of 0-255 values for RGB colour of each CDS
    outline_colours=: list
        3-tuples of 0-255 values for RGB colour of the outline of each CDS
    line_width=: int
        the line width in pixels of the outline of all CDSs
    """

    def __init__(
        self,
        cds_coords=[],
        cds_ids=[],
        cdss=[],
        cds_height=50,
        cds_spacing=50,
        width=800,
        margin=50,
        fill_colours=[],
        outline_colours=[],
        linewidth=2,
    ):
        self.cds_coords = cds_coords
        self.cds_ids = cds_ids
        self.cdss = cdss
        self.cds_height = cds_height
        self.cds_spacing = cds_spacing
        self.width = width
        self.margin = margin
        self.height = (
            (self.cds_height * len(self.cds_coords))
            + (self.cds_spacing * (len(self.cds_coords) - 1))
            + (2 * self.margin)
        )
        self.fill_colours = fill_colours
        self.outline_colours = outline_colours
        self.linewidth = linewidth
        self.img = Image.new("RGB", (self.width, self.height), "white")
        self.canv = ImageDraw.Draw(self.img, mode="RGBA")

    def make_cdss(self):
        """
        Creates a group of CDSs.
        """
        flatten = lambda *n: (
            e
            for a in n
            for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,))
        )
        all_coords = list(flatten(self.cds_coords))
        min_x = min(all_coords)
        max_x = max(all_coords)
        range_x = max_x - min_x
        for i in range(len(self.cds_coords)):
            exon_coords = self.cds_coords[i]
            min_coord = min([coord for coords in exon_coords for coord in coords])
            max_coord = max([coord for coords in exon_coords for coord in coords])
            width_less_margin = self.width - (2 * self.margin)
            x1 = self.margin + (((min_coord - min_x) / range_x) * width_less_margin)
            x2 = (
                self.width
                - self.margin
                - (((max_x - max_coord) / range_x) * width_less_margin)
            )
            y = self.margin + (i * self.cds_height) + (i * self.cds_spacing)
            fill_colour = self.fill_colours[i]
            outline_colour = self.outline_colours[i]
            linewidth = self.linewidth
            if len(self.cds_ids):
                cds_id = self.cds_ids[i]
            else:
                cds_id = None
            cds = CDS(
                self.canv,
                cds_id=cds_id,
                exon_coords=exon_coords,
                x1=x1,
                x2=x2,
                y=y,
                height=self.cds_height,
                fill=fill_colour,
                outline=outline_colour,
                linewidth=linewidth,
            )
            self.cdss.append(cds)
        return self.cdss

    def draw(self, label=False, label_size=16):
        """
        Draws a group of CDSs and returns the resulting
        PIL image object.

            Parameters:
                label (Boolean): whether to add labels to the CDSs (default False)
                label_size (int): the size in points of the label text

            Returns:
                self.img (PIL image object): PIL image with the CDS group drawn

        """
        self.make_cdss()
        max_label_width = 0
        for cds in self.cdss:
            cds.draw()
            if cds.cds_id is not None and label is True:
                label_x = self.width - (0.8 * self.margin)
                label_y = cds.y
                left, top, right, bottom = cds.draw_label(
                    label_x, label_y, label_size=label_size
                )
                label_width = right - left
                max_label_width = max(max_label_width, label_width)
        if label is True:
            remaining_margin = int(self.width - (label_x + max_label_width))
            self.img = add_margin(self.img, right=self.margin - remaining_margin)
        return self.img


def add_margin(img, top=0, right=0, bottom=0, left=0, colour=(255, 255, 255)):
    """
    Takes a PIL image object, adds a margin, and returns a new image.

        Arguments:
            img: a PIL image object created e.g. with Image.new()

        Parameters:
            top (int): number of pixels of margin to add to the top (default 0)
            right (int): number of pixels of margin to add to the right (default 0)
            bottom (int): number of pixels of margin to add to the bottom (default 0)
            left (int): number of pixels of margin to add to the left (default 0)
            colour (3-tuple of ints): 0-255 values for RGB colour of margin (default (255, 255, 255))

        Returns:
            new_img (PIL image object): new PIL image object with margin added
    """
    w, h = img.size
    new_w = w + right + left
    new_h = h + top + bottom
    new_img = Image.new(img.mode, (new_w, new_h), colour)
    new_img.paste(img, (left, top))
    return new_img


if __name__ == "__main__":
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
