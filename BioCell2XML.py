#!/usr/bin/env python3

"""
BioCell2XML.py is a script that converts SIMI BioCell lineage data into the MaMuT-readable XML
file format.

If you use the program for your work please cite the publication:

Pennerstorfer, M., Loose, G. & Wolff, C. (2019). BioCell2XML: A novel tool for converting cell
lineage data from SIMI BioCell to MaMuT (Fiji). Development Genes and Evolution.
https://doi.org/10.1007/s00427-019-00633-9

BioCell2XML reads information from three input files: the SBC and the SBD file of the SIMI
BioCell project that shall be converted, and the H5/XML file of the corresponding image data.

Names of all input files and directories must not contain any space or non-ACSII characters.
Names of the two SIMI BioCell project files have to be identical (beside extensions SBC and SBD).

In order to run BioCell2XML, the image data used in SIMI BioCell have to be converted to the HDF5
container file format used by MaMuT. This will result in an image data file (H5) and in an image
metadata file (XML), the latter is read by BioCell2XML. A detailed description on how this image
pre-processing is done is available in the documentation file BioCell2XML_User_Guide.pdf.

BioCell2XML can be run either in an interactive mode (which will prompt the user step-by-step
through the file and argument input process), or by entering all conversion specifications as
optional arguments directly in the command line.

To run the program in interactive mode enter
$ python3 BioCell2XML.py -p
in the command line. Note, that in this mode, the SIMI BioCell project files (SBC, SBD), the
H5/XML file, and the program file BioCell2XML.py, all have to be located in the same directory.

To run the program in non-interactive mode enter
$ python3 Biocell2XML.py -simi nameofbiocellproject -h5xml nameofh5xmlfile
followed by the optional arguments for specific translation settings (see below).
-simi (followed by the name of the SIMI BioCell project that shall be converted) and -h5xml
(followed by the name of the H5/XML file) are required arguments in this mode. The mode allows
to input these files also from directories which are not the current working directory.

BioCell2XML will output a MaMuT readable XML file with the translated lineage data to the current
working directory.
An additional legend.txt file is written if the user chooses to translate metadata that cannot be
displayed directly in MaMuT (spot shapes, cell fates, user comments), or if the SIMI BioCell lineage
contains alternative branch(es) or mitoses that have only one offspring cell.

List of arguments:
-h      --help          show this help message and exit
-p      --prompt        activate prompt input mode
-simi                   name of SIMI BioCell project (without file extensions .sbc or .sbd)
                        [required argument if prompt input mode is not used]
-h5xml                  name of h5.xml file (without file extension .xml)
                        [required argument if prompt input mode is not used]
-simipath               path of SIMI BioCell project files
-h5xmlpath              path of h5.xml file
-x      --crop_x        crop distance in pixels for x-coordinates
-y      --crop_y        crop distance in pixels for y-coordinates
-z      --crop_z        crop distance in slices for z-coordinates
-m      --mitoses       mode for reconstruction of mitoses (choices: 0, 1, 2, 3):
-m0     --mitoses0      mode 0 (default): no mitosis reconstruction
-m1     --mitoses1      mode 1: offspring spots are added at frame of mitosis
-m2     --mitoses2      mode 2: parent spot is added one frame before mitosis
-m3     --mitoses3      mode 3: offspring spots are added at frame of mitosis and
                                parent spot is added one frame before mitosis
-t      --test          perform test run
-c      --color         deactivate translation of spot colors
-s      --size          deactivate translation of spot sizes
-d      --metadata      deactivate translation of spot shapes, cell fates, and user comments
-i      --interpolate   interpolate between spots
-f      --fraction      fraction for interpolation
-o      --omitmitoses   omit mitoses in interpolation
-out                    name of output file (without file extension .xml)

For a detailed documentation on how use the program see the file BioCell2XML_User_Guide.pdf.
"""

import argparse
import math
import operator
import os
import string
import sys


class Simi:
    """Take in SIMI BioCell project files (SBC, SBD) and store settings and cells."""

    def __init__(self, siminame):
        # Open SIMI BioCell SBC file.
        self.sbc = self.open_file(siminame + '.sbc')
        # Open SIMI BioCell SBD file.
        self.sbd = self.open_file(siminame + '.sbd')

        # Dictionary with all settings from SBC file.
        self.settings = {}
        # Parse SBC data.
        self.parse_sbc()

        # List of all cells that have spot entries from SBD file.
        self.cells = []
        # Parse SBD data.
        self.parse_sbd()

    def open_file(self, filepath):
        """Try to open and return file."""
        try:
            fin = open(filepath, errors='ignore')
            return fin
        except OSError:
            path, file = os.path.split(filepath)
            print("Could not load SIMI BioCell file '{0}' at '{1}'".format(file, path))

    def parse_sbc(self):
        """Parse SBC file and extract all settings."""
        # Skip header of first four lines and read file line by line.
        for line in self.sbc.readlines()[4:]:

            if line.startswith('['):
                title = line.strip('[]\n')
                self.settings[title] = {}

            else:
                key, value = line.strip('\n').split('=')
                self.settings[title][key] = value

    def parse_sbd(self):
        """Parse SBD file and extract all cells that have spots entries."""
        # Skip header of first seven lines.
        sbd = self.sbd.readlines()[7:]

        # Remove possible blank lines in SBD.
        while '\n' in sbd:
            print('blankline in SBD file was removed..')
            sbd.remove('\n')

        # Read file line by line.
        for i in range(len(sbd)-4):

            # Find cell with spot entries.
            if sbd[i].startswith('---') and not sbd[i+5].startswith('---'):

                cell = Cell()

                # Take relevant cell line 1 information and attribute to cell object.
                line = sbd[i+1].split()
                cell.name = line[4]                         # Generic cell name

                # Take relevant cell line 3 information and attribute to cell object.
                line = sbd[i+3].split()
                cell.birth_frame = int(line[0])             # Framenumber of birth of cell (mitosisframe)
                cell.fate = int(line[2])                    # Integer defining cell fate
                cell.size = int(line[3])                    # Integer defining cell size
                cell.shape = int(line[4])                   # Integer defining cell shape
                cell.color = int(line[5])                   # Integer defining cell color
                if len(line) > 6:
                    cell.name = ' '.join(line[6:])          # Custom cell name given by BioCell user, overwrites generic name

                # Take cell line 4 information and attribute to cell object.
                line = sbd[i+4].split()
                cell.nspots = int(line[0])                  # Number of spot entries of cell
                if len(line) > 1:
                    cell.comment = ' '.join(line[1:])       # Biocell user comment for cell

                # Iterate through all spot entries of cell.
                for j in range(5, 5+cell.nspots):

                    spot = Spot()

                    # Take spot line information and attribute to spot object.
                    line = sbd[i+j].split()
                    spot.frame = int(line[0])               # Framenumber of spot
                    spot.x = float(line[1])                 # x-coordinate of spot
                    spot.y = float(line[2])                 # y-coordinate of spot
                    spot.z = float(line[3])                 # z-coordinate of spot
                    spot.size = int(line[4])                # Integer defining spot size
                    spot.shape = int(line[5])               # Integer defining spot shape
                    spot.color = int(line[6])               # Integer defining spot color
                    if len(line) > 7:
                        spot.comment = ' '.join(line[7:])   # BioCell user comment for spot

                    # Attribute cell metadata to spot object if spot metadata are undefined.
                    spot.fate = cell.fate
                    if spot.size == -1:
                        spot.size = cell.size
                    if spot.shape == -1:
                        spot.shape = cell.shape
                    if spot.color == -1:
                        spot.color = cell.color

                    # Combine possible cell and spot comments into spot comment attribute.
                    if spot.comment and not cell.comment:
                        spot.comment = 'Spot: ' + spot.comment
                    elif cell.comment and not spot.comment:
                        spot.comment = 'Cell: ' + cell.comment
                    elif spot.comment and cell.comment:
                        spot.comment = 'Cell: ' + cell.comment + '; Spot: ' + spot.comment

                    cell.spots.append(spot)

                self.cells.append(cell)

class H5xml:
    """Take in H5/XML file and store image and scan data."""

    def __init__(self, h5xmlname):
        # Open H5/XML file.
        self.h5xml_file = self.open_file(h5xmlname + '.xml')

        # Directory with relevant settings from H5/XML file.
        self.settings = {}
        # Parse H5/XML file.
        self.parse_h5xml()

    def open_file(self, filepath):
        """Try to open and return file."""
        try:
            fin = open(filepath)
            return fin
        except OSError:
            path, file = os.path.split(filepath)
            print("Could not load H5/XML file '{0}' at '{1}'".format(file, path))

    def parse_h5xml(self):
        """Parse H5/XML file line by line and extract image and scan data."""
        for line in self.h5xml_file.readlines():

            if line.startswith(' '*8+'<size>'):
                values = line.strip(' <size>/\n').split()
                self.settings['width'] = int(values[0])             # width of image
                self.settings['height'] = int(values[1])            # height of image
                self.settings['nslices'] = int(values[2])           # number of slices in image stack

            if line.startswith(' '*10+'<size>'):
                values = line.strip(' <size>/\n').split()
                self.settings['pixelwidth'] = float(values[0])      # width of voxel
                self.settings['pixelheight'] = float(values[1])     # height of voxel
                self.settings['voxeldepth'] = float(values[2])      # depth of voxel

            if line.startswith(' '*6+'<first>'):
                firstframe = int(line.strip(' <first>/\n'))

            if line.startswith(' '*6+'<last>'):
                lastframe = int(line.strip(' <last>/\n'))
                self.settings['nframes'] = lastframe-firstframe+1   # number of frames

class Cell:
    """Store cell-related information and update frameshift and voxelsize."""

    def __init__(self):
        self.name = ''
        self.birth_frame = None
        self.fate = None
        self.size = None
        self.shape = None
        self.color = None
        self.nspots = None
        self.comment = None

        # List of all spots of cell.
        self.spots = []

    def update_cell(self, simi, h5):
        """Update frameshift of cell, and frameshift and voxelsize of all its spots.

        Input:
        simi -- SIMI BioCell project files data
        h5 -- H5/XML file data
        """

        # Calculate frameshift of track(s) (i.e. start delay from frame 0).
        frameshift = int(simi.settings['DISC']['TIMEOFFSET'])/int(simi.settings['DISC']['SCANTIME'])

        # Perform bankers' rounding of frameshift.
        if (float(frameshift)%1) >= 0.5:
            frameshift = math.ceil(frameshift)
        else:
            frameshift = round(frameshift)

        # Update birth of cell by frameshift.
        self.birth_frame -= frameshift

        # Calculate timeinterval (i.e. time in seconds).
        timeinterval = int(simi.settings['DISC']['SCANTIME'])/10

        # Update each spot in cell.
        for spot in self.spots:

            # Update temporal attributes by frameshift.
            spot.frame -= frameshift
            spot.t -= frameshift*timeinterval

            # Correct voxelsize dimension.
            spot.x *= h5.settings['width']/int(simi.settings['CALIBRATION']['WIDTH'])
            spot.y *= h5.settings['width']/int(simi.settings['CALIBRATION']['WIDTH'])
            spot.z *= h5.settings['voxeldepth']/h5.settings['pixelwidth']

class Spot:
    """Store spot-related information.

    The class provides functions to add spot attributes, convert color and size from BioCell's
    format to the format used in MaMuT, crop spatial coordinates, interpolate between spots,
    and calculate the intermediate position between spots.
    """

    def __init__(self):
        self.ID = 0
        self.name = ''
        self.frame = None
        self.x = None
        self.y = None
        self.z = None
        self.t = None
        self.fate = None
        self.size = None
        self.shape = None
        self.color = None
        self.comment = None

    def add_spot_attributes(self, spotID, cell, simi, h5, args):
        """Add and convert spot attributes to format used in XML file.

        Input:
        spotID -- next assignable spot ID (integer)
        cell -- cell object to which the spot belongs
        simi -- SIMI BioCell project files data
        h5 -- H5/XML file data
        args -- program arguments

        Return:
        spotID -- next assignable spot ID (integer)
        """

        # Add spotID and name attributes to spot object.
        self.ID = spotID
        self.name = cell.name + '_ID' + str(self.ID)
        spotID += 1

        # Generate t attribute of spot object.
        timeinterval = int(simi.settings['DISC']['SCANTIME'])/10
        self.t = float(self.frame)*timeinterval

        # Convert spot size (if defined in BioCell), or set size to MaMuT's default value.
        if args.size and self.size != -1:
            self.convert_size(self.size, h5)
        else:
            self.size = 10.0

        # Convert spot color (if defined in BioCell).
        if args.color and self.color != -1:
            self.convert_color(self.color)
        else:
            self.color = None

        return spotID

    def convert_size(self, size, h5):
        """Convert spot size from BioCell's format to resolution-dependent format used by MaMuT."""
        self.size = ((0.01*size + 0.06) * h5.settings['width'])/2

    def convert_color(self, ole):
        """Convert spot color from BioCell's OLE format to WIN format used by MaMuT."""
        # Calculate value for each color channel (8-bit).
        blue = (ole/65536) % 256
        green = (ole/256) % 256
        red = ole % 256

        # Convert RGB to WIN format (decimal).
        win = int(red*65536 + green*256 + blue - 16777216)

        # Deal with exceptional OLE value 65535 (which represents yellow (#FFFF00),
        # but results in a WIN value of 0).
        if ole == 65535:
            self.color = -256
        else:
            self.color = win

    def crop_coordinates(self, args):
        """Crop spot coordinates in x, y, and z dimension."""
        self.x -= args.crop_x
        self.y -= args.crop_y
        self.z -= args.crop_z

    def interpolate(self, other, fraction=1, nframes=None, spotID=None, targetframe=None):
        """Calculate new spot(s) by linear interpolation between two input spots.

        The function generates either a series of new spot objects (at every multiple framenumber
        of fraction) or, if a targetframe is given, a single new spot object at that targetframe.

        Input:
        other -- spot object
        fraction -- frame distance between interpolated spots (integer)

        Optional arguments:
        nframes -- number of frames (integer)
        spotID -- next assignable spot ID (integer)
        targetframe -- framenumber at which a spot will be interpolated (integer)

        Return:
        new_spot -- single spot object
        interpolated -- list of spot objects
        spotID -- next assignable spot ID (integer)
        """

        assert fraction > 0

        # Calculate step distances for interpolation.
        frame_difference = other.frame - self.frame
        step_t = (other.t-self.t) / frame_difference
        step_x = (other.x-self.x) / frame_difference
        step_y = (other.y-self.y) / frame_difference
        step_z = (other.z-self.z) / frame_difference
        step_size = (other.size-self.size) / frame_difference

        # Set start if targetframe is given.
        if targetframe is not None:
            start = 1

        # Set start to first multiple of fraction that comes after framenumber of first spot.
        else:
            if self.frame/fraction == int:
                start = fraction
            else:
                for x in range(fraction, nframes, fraction):
                    if x > self.frame:
                        start = x-self.frame
                        break

            # List for series of interpolated spots.
            interpolated = []

        # Iterate through frames between the two spots, beginning at start variable and
        # advancing with step distance fraction.
        for i in range(start, frame_difference, fraction):

            # Interpolate only one spot at targetframe.
            if targetframe is not None:
                i = targetframe - self.frame

            new_spot = Spot()

            # Interpolate frame, t, x, y, z values and attribute to new spot object.
            new_spot.frame = self.frame + i
            new_spot.t = self.t + step_t*i
            new_spot.x = self.x + step_x*i
            new_spot.y = self.y + step_y*i
            new_spot.z = self.z + step_z*i

            # Return that one interpolated spot.
            if targetframe is not None:
                return new_spot

            # Add spotID and name attributes to spot.
            new_spot.ID = spotID
            spotID += 1

            # In case of mitosis take name from offspring cell.
            if self.name.split('ID')[0] != other.name.split('ID')[0]:
                new_spot.name = other.name.split('ID')[0] + 'ID' + str(new_spot.ID)
            else:
                new_spot.name = self.name.split('ID')[0] + 'ID' + str(new_spot.ID)

            # Interpolate and attribute spot size.
            new_spot.size = self.size + step_size*i

            # Set color of interpolated spot.
            if i <= frame_difference/2:
                new_spot.color = self.color
            else:
                new_spot.color = other.color

            interpolated.append(new_spot)

        return interpolated, spotID

    def intermediate(self, others):
        """Calculate and return a spot at the intermediate position between the input spots.

        Input:
        others -- list of spot objects
        """

        new_spot = Spot()
        new_spot.frame = self.frame
        new_spot.t = self.t

        # For x, y, z sum-up values of all input spots.
        x = self.x
        y = self.y
        z = self.z

        for spot in others:
            x += spot.x
            y += spot.y
            z += spot.z

        # Calculate center x, y, z coordinates and attribute to new spot object.
        new_spot.x = x / (1+len(others))
        new_spot.y = y / (1+len(others))
        new_spot.z = z / (1+len(others))

        return new_spot

class Mitosis():
    """Store mitosis-related information."""

    def __init__(self, parentspot, splitframe):
        # Original parent spot of mitosis in BioCell (i.e. last spot entry of a cell before its
        # division in BioCell lineage data).
        self.parent = parentspot

        # Framenumber of birth of offspring cell (mitosisframe) (Integer).
        self.splitframe = splitframe

        # List of original offspring spots in BioCell (i.e. first spot entry in offspring cell
        # in BioCell lineage data).
        self.offsprings = []

        # Dictionary mapping original parent spot to new parent spot (coming from function
        # calculate_mitosis_spots in class Track).
        self.new_parent = {}

        # Dictionary mapping each original offspring spot to new offspring spot (coming from
        # function calculate_mitosis_spots in class Track).
        self.new_offsprings = {}

class Edge:
    """An edge is a connection between two spot objects (sourcespot and targetspot)."""

    def __init__(self, sourcespot, targetspot):
        self.source = sourcespot
        self.target = targetspot

    def duration(self):
        """Return frame distance between source and target."""
        return float(self.target.frame-self.source.frame)

    def displacement(self):
        """Return displacement between source and target."""
        return math.sqrt((self.source.x-self.target.x)**2 + (self.source.y-self.target.y)**2 + (self.source.z-self.target.z)**2)

    def velocity(self):
        """Return velocity of edge."""
        return self.displacement()/self.duration()

class Track():
    """A track is (the MaMuT terminology of) a cell lineage.

    The class reconstructs the track out of the BioCell SBD data, generates and stores track data,
    and provides functions to calculate new mitosis spots, replace original mitosis edges by edges
    between the new mitosis spots, and interpolate new spots in the track.
    """

    def __init__(self, trackID, cells):
        """Input:
        trackID -- consecutive tracknumber (Integer)
        cells -- list of cells from BioCell SBD file
        """

        self.ID = trackID
        self.name = 'Track_' + str(self.ID)

        # List of all spots in track.
        self.spots = []
        # List of all edges in track.
        self.edges = []
        # Number of splits in track.
        self.nsplits = 0

        # Number of spots in track.
        self.nspots = 0
        # Number of gaps between the spots.
        self.ngaps = 0

        # Framenumber of first spot in track.
        self.start = None
        # Framenumber of last spot in track.
        self.stop = None
        # Frame distance between first and last spot.
        self.duration = 0
        # Displacement between first and last spot.
        self.displacement = 0
        # Frame distance of edge that has the longest distance between source and target spot.
        self.longestgap = 0

        # Mapping of parent spot to mitosis object.
        self.mitoses = {}
        # List of all edges that belong to a mitosis.
        self.mitosesedges = []

        # Index of last cell of track in cells list.
        self.end = 0

        # Reconstruct track.
        self.reconstruct_track(cells)

    def reconstruct_track(self, cells):
        """Reconstruct cell genealogy and find end of a track.

        The function iterates through the list of cells from BioCell SBD data, finds mitosis events
        and alternative branches, adds spots and edges to the track data, and finds the end of a
        track in the cells list.

        Input:
        cells -- list of cells from BioCell SBD file
        """

        # Add edges for first cell in cells.
        self.add_edges(cells[0].spots)

        # End of track (in case cells contains only one cell).
        if len(cells) == 1:
            self.end = 1
            return

        # Iterate through cells.
        for i in range(1, len(cells)):

            offspring_spot = cells[i].spots[0]

            # Framenumber of birth of cell (mitosisframe).
            birth_cell = cells[i].birth_frame
            # Framenumber of birth of previous cell in cells.
            birth_previous_cell = cells[i-1].birth_frame

            if birth_cell > birth_previous_cell:

                # End of track.
                if birth_cell < cells[i-1].spots[-1].frame:
                    self.end = i
                    return

                # Mitosis event.
                parent_spot = cells[i-1].spots[-1]
                self.add_mitosis(parent_spot, offspring_spot, birth_cell)

            elif birth_cell == birth_previous_cell:

                # End of track.
                if i-2 < 0:
                    self.end = i
                    return

                # Mitosis event.
                if birth_cell > cells[i-2].spots[-1].frame:
                    parent_spot = cells[i-2].spots[-1]
                    self.add_mitosis(parent_spot, offspring_spot, birth_cell)

                # Find alternative branch.
                else:
                    # Iterate backwards through cells.
                    for j in range(i-3, -1, -1):
                        if birth_cell == cells[j].birth_frame and birth_cell > cells[j-1].spots[-1].frame:
                            parent_spot = cells[j-1].spots[-1]
                            self.add_mitosis(parent_spot, offspring_spot, birth_cell)

            elif birth_cell < birth_previous_cell:

                # Iterate backwards through cells.
                for j in range(i-1, -1, -1):

                    # End of track.
                    if j == 0:
                        self.end = i
                        return

                    # Mitosis event or alternative branch.
                    if birth_cell == cells[j].birth_frame and birth_cell > cells[j-1].spots[-1].frame:
                        parent_spot = cells[j-1].spots[-1]
                        self.add_mitosis(parent_spot, offspring_spot, birth_cell)
                        break

            # Add edges for offspring cell.
            self.add_edges(cells[i].spots)

        # End of track.
        self.end = i+1
        return

    def add_edges(self, spots, index=None):
        """Add spots and edges between them to respective lists of track.

        If index is given insert edge at this position in edges list.

        Input:
        spots -- list of spot objects

        Optional arguments:
        index -- position for insertion (integer)
        """

        for i in range(len(spots)):
            spot = spots[i]

            if spot not in self.spots:
                self.spots.append(spot)

            if i == len(spots)-1:
                break

            edge = Edge(spot, spots[i+1])

            if index is not None:
                self.edges.insert(index, edge)
                index += 1
                self.mitosesedges.append(edge)
            else:
                self.edges.append(edge)

    def add_mitosis(self, parent, offspring, splitframe):
        """Instantiate and add mitosis object to mitoses list, increase number of splits, and
        add mitosis edges to lists edges and mitosesedges.

        Input:
        parent -- parent spot object of mitosis event
        offspring -- offspring spot object of mitosis event
        splitframe -- framenumber of birth of offspring cell (mitosisframe) (integer)
        """

        if parent not in self.mitoses:
            mitosis = Mitosis(parent, splitframe)
            mitosis.offsprings.append(offspring)
            self.mitoses[parent] = mitosis

        else:
            self.mitoses[parent].offsprings.append(offspring)

        self.nsplits += 1

        edge = Edge(parent, offspring)
        self.edges.append(edge)
        self.mitosesedges.append(edge)

    def calculate_mitosis_spots(self, spotID, args):
        """Reconstruct new mitosis (parent and/or offspring) spots.

        The new parent spot is reconstructed by linear interpolation between the original parent
        spot and each of the original offspring spots (from BioCell lineage data) onto the frame
        one framenumber before the birth of the offspring cell (mitosisframe-1), and subsequent
        calculation of the intermediate position between these interpolated spots.
        New offspring spots are reconstructed by interpolation between the original parent spot
        and the original offspring spot onto the birth frame of the offspring cell (mitosisframe).

        If the program argument --test is True the color of new mitosis spots is set to green,
        or, if the new mitosis spot coincides with the original mitosis spot, to magenta.
        If --test is False spot color and size are adopted from the original mitosis spot.

        Input:
        spotID -- next assignable spot ID (integer)
        args -- program arguments

        Return:
        spotID -- next assignable spot ID (integer)
        """

        # Iterate through mitoses dictionary.
        for parent, mitosis in self.mitoses.items():

            # Reconstruct new parent spot.
            if args.mitoses == 2 or args.mitoses == 3:

                # Original parent spot lies on mitosisframe-1.
                if parent.frame == mitosis.splitframe-1:

                    if args.test:
                        # Set spot color to magenta.
                        parent.color = -65281

                else:
                    new_spots = []

                    # For each offspring interpolate a new spot onto mitosisframe-1.
                    for offspring in mitosis.offsprings:
                        new_spot = parent.interpolate(offspring, targetframe=mitosis.splitframe-1)
                        new_spots.append(new_spot)

                    # Generate new parent spot as geometic center between these new spots.
                    new_parent = new_spots[0].intermediate(new_spots[1:])

                    # Attribute ID and name to new parent spot.
                    new_parent.ID = spotID
                    new_parent.name = parent.name.split('ID')[0] + 'ID' + str(new_parent.ID) + '_mitosis'
                    spotID += 1

                    if args.test:
                        # Set spot color to green and spot size to MaMuT's default value.
                        new_parent.color = -16711936
                        new_parent.size = 10.0

                    else:
                        # Adopt spot color and spot size from original parent spot.
                        new_parent.color = parent.color
                        new_parent.size = parent.size

                    # Add new parent spot to mitosis object.
                    self.mitoses[parent].new_parent[parent] = new_parent

            # Reconstruct new offspring spots.
            if args.mitoses == 1 or args.mitoses == 3:

                # For each offspring of the mitosis.
                for offspring in mitosis.offsprings:

                    # Offspring spot lies on mitosisframe.
                    if offspring.frame == mitosis.splitframe:

                        if args.test:
                            # Set spot color to magenta.
                            offspring.color = -65281

                    # Interpolate new offspring spot onto mitosisframe.
                    else:
                        new_offspring = parent.interpolate(offspring, targetframe=mitosis.splitframe)

                        # Attribute ID and name to new offspring spot.
                        new_offspring.ID = spotID
                        new_offspring.name = offspring.name.split('ID')[0] + 'ID' + str(new_offspring.ID) + '_mitosis'
                        spotID += 1

                        if args.test:
                            # Set spot color to green and spot size to MaMuT's default value.
                            new_offspring.color = -16711936
                            new_offspring.size = 10.0

                        else:
                            # Adopt spot color and spot size from original offspring spot.
                            new_offspring.color = offspring.color
                            new_offspring.size = offspring.size

                        # Add new offspring spot to mitosis object.
                        self.mitoses[parent].new_offsprings[offspring] = new_offspring

        return spotID

    def replace_mitoses_edges(self, args):
        """Replace original mitoses edges by new mitoses edges.

        The function deletes edges between original parent and offspring spots (from BioCell
        lineage data) from the track, and adds new mitosis spots (from function
        calculate_mitosis_spots) and the corresponding mitosis edges to the track.

        Depending on the chosen program argument --mitoses either only new offsping spots, only
        new parent spots, or both, new offspring spots and parent spots, are added.

        Input:
        args -- program arguments
        """

        # Iterate through edges list.
        for edge in self.edges:

            if edge.source in self.mitoses:

                # Get index of original edge in edges list.
                index = self.edges.index(edge)

                # Reconstruction mode 1: Add new offspring spots at mitosisframe.
                if args.mitoses == 1:

                    if edge.target in self.mitoses[edge.source].new_offsprings:

                        new_target = self.mitoses[edge.source].new_offsprings[edge.target]

                        # Remove original edge.
                        self.edges.remove(edge)
                        self.mitosesedges.remove(edge)

                        # Add new mitosis spot and new mitosis edges.
                        self.add_edges([edge.source, new_target, edge.target], index=index)

                # Reconstruction mode 2: Add new parent spots at mitosisframe-1.
                if args.mitoses == 2:

                    if edge.source in self.mitoses[edge.source].new_parent:

                        new_source = self.mitoses[edge.source].new_parent[edge.source]

                        # Remove original edge.
                        self.edges.remove(edge)
                        self.mitosesedges.remove(edge)

                        # Add new mitosis spot and new mitosis edges.
                        # If first split event in mitosis.
                        if edge.target == self.mitoses[edge.source].offsprings[0]:
                            self.add_edges([edge.source, new_source, edge.target], index=index)

                        # If further split events in mitosis.
                        else:
                            self.add_edges([new_source, edge.target], index=index)

                # Reconstruction mode 3: Add new parent spots at mitosisframe-1 and
                # new offspring spots at mitosisframe.
                if args.mitoses == 3:

                    if edge.source in self.mitoses[edge.source].new_parent and edge.target in self.mitoses[edge.source].new_offsprings:

                        new_source = self.mitoses[edge.source].new_parent[edge.source]
                        new_target = self.mitoses[edge.source].new_offsprings[edge.target]

                        # Remove original edge.
                        self.edges.remove(edge)
                        self.mitosesedges.remove(edge)

                        # Add new mitosis spots and new mitosis edges.
                        # If first split event in mitosis.
                        if edge.target == self.mitoses[edge.source].offsprings[0]:
                            self.add_edges([edge.source, new_source, new_target, edge.target], index=index)

                        # If further split events in mitosis.
                        else:
                            self.add_edges([new_source, new_target, edge.target], index=index)

                    elif edge.source in self.mitoses[edge.source].new_parent and edge.target not in self.mitoses[edge.source].new_offsprings:

                        new_source = self.mitoses[edge.source].new_parent[edge.source]

                        # Remove original edge.
                        self.edges.remove(edge)
                        self.mitosesedges.remove(edge)

                        # Add new mitosis spot and new mitosis edges.
                        # If first split event in mitosis.
                        if edge.target == self.mitoses[edge.source].offsprings[0]:
                            self.add_edges([edge.source, new_source, edge.target], index=index)

                        # If further split events in mitosis.
                        else:
                            self.add_edges([new_source, edge.target], index=index)

                    elif edge.target  in self.mitoses[edge.source].new_offsprings and edge.source not in self.mitoses[edge.source].new_parent:

                        new_target = self.mitoses[edge.source].new_offsprings[edge.target]

                        # Remove original edge.
                        self.edges.remove(edge)
                        self.mitosesedges.remove(edge)

                        # Add new mitosis spot and new mitosis edges.
                        self.add_edges([edge.source, new_target, edge.target], index=index)

    def interpolate_track(self, spotID, nframes, total, runs, args):
        """Add spots between existing spots by linear interpolation.

        The function adds new spots and edges between them to the track, and subsequently deletes
        all previous edges from the track's edges list. A console progress bar for this process
        is updated.

        If the program argument --omitmitoses is True no interpolation is done between parent spot
        and offspring spots of the original BioCell lineage.

        Input:
        spotID -- next assignable spot ID (integer)
        nframes -- number of frames (integer)
        total -- total number of edges in all tracks (integer)
        runs -- number of processed edges (integer)
        args -- program arguments

        Return:
        spotID -- next assignable spot ID (integer)
        runs -- number of processed edges (integer)
        """

        end = len(self.edges)

        # Iterate through edges list.
        for edge in self.edges[:end]:

            # Omit edges that belong to a mitosis.
            if args.omitmitoses and edge in self.mitosesedges:
                self.edges.append(edge)

            # Interpolate spots and add new spots and edges to track.
            else:
                interpolated, spotID = edge.source.interpolate(edge.target, fraction=args.fraction, nframes=nframes, spotID=spotID)
                interpolated.insert(0, edge.source)
                interpolated.append(edge.target)
                self.add_edges(interpolated)

            # Update progress bar.
            runs += 1
            update_progress(runs/total)

        # Delete all previous edges from track.
        del self.edges[:end]

        return spotID, runs

    def get_track_data(self):
        """Get other track data."""

        self.nspots = len(self.spots)
        self.ngaps = len(self.edges)

        # First spot in track.
        firstspot = self.spots[0]
        # Last spot in track (i.e. spot with highest framenumber).
        lastspot = max(self.spots, key=operator.attrgetter('frame'))

        self.start = firstspot.frame
        self.stop = lastspot.frame

        self.duration = Edge(firstspot, lastspot).duration()
        self.displacement = Edge(firstspot, lastspot).displacement()

        # Find longest gap in track.
        for edge in self.edges:
            gap = edge.duration()
            if gap > self.longestgap:
                self.longestgap = gap


def biocell2xml():
    """Convert SIMI BioCell lineage data to MaMuT-readable XML file format."""

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Converts SIMI BioCell cell lineage data to XML format.\n(ATTENTION: All file and directory names must not contain any space or non-ASCII characters)', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-p', '--prompt', action='store_true', default=False, help='activate prompt input mode')
    parser.add_argument('-simi', metavar='', help='name of SIMI BioCell project (without file extensions .sbc or .sbd)')
    parser.add_argument('-h5xml', metavar='', help='name of h5.xml file (without file extension .xml)')
    parser.add_argument('-simipath', metavar='', help='path of SIMI BioCell project files')
    parser.add_argument('-h5xmlpath', metavar='', help='path of h5.xml file')
    parser.add_argument('-x', '--crop_x', metavar='', type=int, default=0, help='crop distance in pixels for x-coordinates')
    parser.add_argument('-y', '--crop_y', metavar='', type=int, default=0, help='crop distance in pixels for y-coordinates')
    parser.add_argument('-z', '--crop_z', metavar='', type=int, default=0, help='crop distance in slices for z-coordinates')
    parser.add_argument('-m', '--mitoses', metavar='', type=int, default=0, choices=[0, 1, 2, 3], help='mode for reconstruction of mitoses (choices: 0, 1, 2, 3)')
    parser.add_argument('-t', '--test', action='store_true', default=False, help='perform test run')
    parser.add_argument('-c', '--color', action='store_false', default=True, help='deactivate translation of spot colors')
    parser.add_argument('-s', '--size', action='store_false', default=True, help='deactivate translation of spot sizes')
    parser.add_argument('-d', '--metadata', action='store_false', default=True, help='deactivate translation of spot shapes, cell fates, and user comments')
    parser.add_argument('-i', '--interpolate', action='store_true', default=False, help='interpolate between spots')
    parser.add_argument('-f', '--fraction', metavar='', type=int, default=1, help='fraction for interpolation')
    parser.add_argument('-o', '--omitmitoses', action='store_true', default=False, help='omit mitoses in interpolation')
    parser.add_argument('-out', metavar='', help='name of output file (without file extension .xml)')

    args = parser.parse_args()

    # Call function to input program arguments via prompting the user.
    if args.prompt:
        args = prompt_input(args)

    # Take in data from BioCell SBC and SBD files.
    if args.simipath:
        simi = Simi(os.path.join(args.simipath, args.simi))
    else:
        simi = Simi(os.path.join(os.getcwd(), args.simi))

    # Take in data from H5/XML file.
    if args.h5xmlpath:
        h5 = H5xml(os.path.join(args.h5xmlpath, args.h5xml))
    else:
        h5 = H5xml(os.path.join(os.getcwd(), args.h5xml))

    # For test run (i.e. if argument --test is True) override program arguments.
    if args.test:
        args.color = False
        args.size = False
        args.metadata = False
        args.interpolate = False

    # Transform and update spot attributes to format used in XML file.
    spotID = 1
    for cell in simi.cells:
        for spot in cell.spots:
            spotID = spot.add_spot_attributes(spotID, cell, simi, h5, args)
        cell.update_cell(simi, h5)

    # Crop spot coordinates.
    for cell in simi.cells:
        for spot in cell.spots:
            spot.crop_coordinates(args)

    # Reconstruct tracks and store them in the list tracks.
    tracks = []
    trackID = 1
    # Iterate through cells from SBD file until all tracks are reconstructed.
    i = 0
    while len(simi.cells[i:]) > 0:
        track = Track(trackID, simi.cells[i:])
        tracks.append(track)
        trackID += 1
        i += track.end

    # Generate new mitosis spots and replace mitosis edges.
    if args.mitoses != 0:
        for track in tracks:
            spotID = track.calculate_mitosis_spots(spotID, args)
            track.replace_mitoses_edges(args)

    # Interpolate between spots.
    if args.interpolate:

        # Data for progress bar.
        total = 0
        for track in tracks:
            total += len(track.edges)
        runs = 0

        # Call interpolation function.
        for track in tracks:
            spotID, runs = track.interpolate_track(spotID, int(h5.settings['nframes']), total, runs, args)

    # Output XML file to current working directory.
    xml_filename, flag_shape, flag_fate, comments = output_xml_file(tracks, simi, h5, args)

    # Output legend.txt file to current working directory.
    legend_filename = output_legend_file(tracks, simi, xml_filename, flag_shape, flag_fate, comments, args)

    print('Transformation completed:')
    print("File '{0}' was written to '{1}'".format(xml_filename, os.getcwd()))
    if legend_filename is not None:
        print("File '{0}' was written to '{1}'".format(legend_filename, os.getcwd()))

def prompt_input(args):
    """Prompt user to input program arguments step-by-step and return them as args."""

    print('ATTENTION: BioCell2XML.py, the SIMI BioCell SBC and SBD files, and the H5.XML file must be located in the same directory\nATTENTION: All file and directory names must not contain any space or non-ASCII characters')

    args.simi = input('Enter name of SIMI BioCell project (without file extensions .sbc or .sbd): ')
    args.h5xml = input('Enter name of h5.xml file (without file extension .xml): ')

    crop = input('Did you crop the original image stack? (y/n): ')
    if crop == 'y':
        args.crop_x = int(input('Enter crop distance in pixels for x-coordinates (integer value; will be subtracted from left image border): '))
        args.crop_y = int(input('Enter crop distance in pixels for y-coordinates (integer value; will be subtracted from upper image border): '))
        args.crop_z = int(input('Enter crop distance in slices for z-coordinates (integer value; will be subtracted from lower border of image stack): '))

    args.mitoses = int(input('Choose mode for reconstruction of mitoses:\n0: no mitosis reconstruction\n1: additional offspring spots at frame of mitosis\n2: additional parent spot one frame before mitosis\n3: additional offspring spots at frame of mitosis and additional parent spot one frame before mitosis\n(0, 1, 2, or 3): '))

    test = input('Perform test run (this option will highlight all mitoses spots)? (y/n): ')
    if test == 'y':
        args.test = True
        return args

    metadata = input('Translate SIMI BioCell spot metadata to XML file? (y/n): ')
    if metadata == 'n':
        args.color, args.size, args.metadata = False, False, False
    elif metadata == 'y':
        color = input('Translate spot colors? (y/n): ')
        if color == 'n':
            args.color = False
        size = input('Translate spot sizes? (y/n): ')
        if size == 'n':
            args.size = False
        metadata = input('Translate spot shapes, cell fates, and user comments? (y/n): ')
        if metadata == 'n':
            args.metadata = False

    interpolate = input('Interpolate between spots (i.e. adding mean coordinate spots between existing spots)? (y/n): ')
    if interpolate == 'y':
        args.interpolate = True
        args.fraction = int(input('Enter fraction (i.e. distance between spots) for interpolation (integer value): '))
        omitmitoses = input('Omit mitoses in interpolation? (y/n): ')
        if omitmitoses == 'y':
            args.omitmitoses = True

    args.out = input('Enter name of output file (optional; without file extension .xml): ')

    return args

def update_progress(progress):
    """Display and update a console progress bar.

    The function is modified after Brian Khuu from:
    https://stackoverflow.com/questions/3160699/python-progress-bar
    """

    bar_length = 20
    status = ''

    if progress == 1:
        status = ' \r\n'

    block = int(round(bar_length*progress))
    text = '\rInterpolation: [{0}] {1}% {2}'.format('#'*block + '-'*(bar_length-block), round(progress*100, 4), status)
    sys.stdout.write(text)
    sys.stdout.flush()

def output_xml_file(tracks, simi, h5, args):
    """Write XML file to current working directory and return filename and metadata flags."""

    # Define filename.
    if args.out:
        xml_filename = args.out + '.xml'
    elif args.test:
        xml_filename = args.simi + '_BioCell2XML_test.xml'
    else:
        xml_filename = args.simi + '_BioCell2XML.xml'

    fout = open(xml_filename, 'w')

    # Write XML header information.
    fout.write(BEGIN_TEMPLATE_1)

    if args.metadata:
        fout.write(BEGIN_METADATA_TEMPLATE)

    fout.write(BEGIN_TEMPLATE_2)

    # Generate a list with all spots of all tracks.
    spots = []
    for track in tracks:
        spots.extend(track.spots)
    # Sort this list by increasing framenumbers.
    spots.sort(key=operator.attrgetter('frame'))

    # Write Spots section.
    fout.write(ALLSPOTS_BEGIN_TEMPLATE.format(nspots=len(spots)))

    filt = None
    flag_shape = False
    flag_fate = False
    comments = {}

    for i in range(len(spots)):
        spot = spots[i]

        # Condition for first spot with a new framenumber.
        if spot.frame != filt:
            fout.write(INFRAME_BEGIN_TEMPLATE.format(frame=spot.frame))
        filt = spot.frame

        # Write spot line.
        fout.write(SPOT_TEMPLATE.format(id=spot.ID, name=spot.name, frame=spot.frame, t=spot.t, x=spot.x, y=spot.y, z=spot.z, radius=spot.size))

        # Write spot color information.
        if spot.color:
            fout.write(SPOT_COLOR_TEMPLATE.format(color=spot.color))

        # Write spot metadata information.
        if args.metadata:

            if spot.shape and spot.shape != -1:
                fout.write(SPOT_SHAPE_TEMPLATE.format(shape=spot.shape))
                flag_shape = True

            if spot.fate and spot.fate != -1:
                fout.write(SPOT_FATE_TEMPLATE.format(fate=spot.fate))
                flag_fate = True

            if spot.comment:
                fout.write(SPOT_COMMENT_TEMPLATE.format(id=spot.ID))
                comments[spot.ID] = spot.comment

        fout.write('/>\n')

        # Condition for last spot with that framenumber.
        try:
            nextspot = spots[i+1]
            if nextspot.frame != filt:
                fout.write(INFRAME_END_TEMPLATE)
        except IndexError:
            fout.write(INFRAME_END_TEMPLATE)

    fout.write(ALLSPOTS_END_TEMPLATE)

    # Write Tracks section.
    fout.write(ALLTRACKS_BEGIN_TEMPLATE)

    for track in tracks:

        track.get_track_data()

        # Write track line.
        fout.write(TRACK_BEGIN_TEMPLATE.format(name=track.name, id=track.ID, duration=float(track.duration), start=float(track.start), stop=float(track.stop), displacement=track.displacement, nspots=track.nspots, ngaps=track.ngaps, longestgap=track.longestgap, nsplits=track.nsplits))

        # Write edge lines.
        for edge in track.edges:
            fout.write(EDGE_TEMPLATE.format(source=edge.source.ID, target=edge.target.ID, velocity=edge.velocity(), displacement=edge.displacement()))

            # Write edge color information.
            if edge.target.color:
                fout.write(EDGE_COLOR_TEMPLATE.format(color=edge.target.color))

            fout.write('/>\n')

        fout.write(TRACK_END_TEMPLATE)

    fout.write(ALLTRACKS_END_TEMPLATE)

    # Write FilteredTracks section.
    fout.write(FILTEREDTRACKS_BEGIN_TEMPLATE)
    for track in tracks:
        fout.write(FILTEREDTRACKS_TEMPLATE.format(id=track.ID))
    fout.write(FILTEREDTRACKS_END_TEMPLATE)

    # Write SIMIBioCellLegends section.
    if flag_shape or flag_fate or comments:
        fout.write(BIOCELLLEGENDS_BEGIN_TEMPLATE)

        if args.metadata:

            # Write spot shapes legend.
            if flag_shape:
                fout.write(SPOT_SHAPES_TEMPLATE)

            # Write cell fates legend.
            if flag_fate:
                for key, value in simi.settings['FATE'].items():
                    index = key.lstrip('FATE')
                    fate = value.rstrip('!').replace(' ', '_').split(',')[1]
                    fout.write(CELL_FATES_TEMPLATE.format(fate=fate, index=index))

            # Write BioCell user comments legend.
            for spotID, comment in sorted(comments.items()):
                # Omit non-ACSII characters and replace characters " and < in comments, as
                # these characters cause problems when reading the XML file with MaMuT.
                for character in comment:
                    if character not in string.printable:
                        comment = comment.replace(character, '')
                    if character == '"':
                        comment = comment.replace('"', "'")
                    if character == '<':
                        comment = comment.replace('<', "'smaller than'")

                fout.write(USER_COMMENT_TEMPLATE.format(id=spotID, comment=comment))

        fout.write(BIOCELLLEGENDS_END_TEMPLATE)

    # Write ImageData.
    if args.h5xmlpath:
        xml_folder = args.h5xmlpath
    else:
        xml_folder = os.getcwd()

    fout.write(IMAGE_DATA_TEMPLATE.format(file=args.h5xml+'.xml', dir=xml_folder, width=h5.settings['width'], height=h5.settings['height'], nslices=h5.settings['nslices'], nframes=h5.settings['nframes'], pixelwidth=h5.settings['pixelwidth'], pixelheight=h5.settings['pixelheight'], voxeldepth=h5.settings['voxeldepth'], timeinterval=int(simi.settings['DISC']['SCANTIME'])/10))

    # Write XML ending information.
    fout.write(END_TEMPLATE)

    fout.close()

    return xml_filename, flag_shape, flag_fate, comments

def output_legend_file(tracks, simi, xml_filename, flag_shape, flag_fate, comments, args):
    """Write legend.txt file to current working directory and return filename."""

    # Find polytomies and mitoses with only one offspring cell.
    single_offspring_mitoses = []
    multiple_offspring_mitoses = []

    for track in tracks:
        for parent, mitosis in track.mitoses.items():
            if len(mitosis.offsprings) == 1:
                single_offspring_mitoses.append(parent.name.split('_ID')[0])
            if len(mitosis.offsprings) > 2:
                multiple_offspring_mitoses.append(parent.name.split('_ID')[0])

    if single_offspring_mitoses or multiple_offspring_mitoses or args.metadata and flag_shape or flag_fate or comments:

        # Define filename.
        legend_filename = xml_filename.replace('.xml', '_legend.txt')

        # Output file.
        fout = open(legend_filename, 'w')

        fout.write(TITLE_TEMPLATE.format(file=xml_filename))

        # Write odd mitoses section.
        if single_offspring_mitoses or multiple_offspring_mitoses:
            fout.write(ODD_MITOSES_HEAD)

            if single_offspring_mitoses:
                fout.write(ONE_OFFSPRING_HEAD)
                for cellname in single_offspring_mitoses:
                    fout.write(ODD_MITOSES_TEMPLATE.format(name=cellname))
                fout.write('\n')

            if multiple_offspring_mitoses:
                fout.write(MULTIPLE_OFFSPRING_HEAD)
                for cellname in multiple_offspring_mitoses:
                    fout.write(ODD_MITOSES_TEMPLATE.format(name=cellname))
                fout.write('\n')

            fout.write('\n')

        # Write spot shapes legend.
        if flag_shape:
            fout.write(SHAPES_TEMPLATE)

        # Write cell fates legend.
        if flag_fate:
            fout.write(FATES_HEAD_TEMPLATE)
            for key, value in simi.settings['FATE'].items():
                index = key.lstrip('FATE')
                fate = value.rstrip('!').split(',')[1]
                fout.write(FATES_TEMPLATE.format(fate=fate, index=index))
            fout.write('\n')

        # Write user comments legend.
        if comments:
            fout.write(COMMENT_HEAD_TEMPLATE)
            for spotID, comment in sorted(comments.items()):
                fout.write(COMMENT_TEMPLATE.format(id=spotID, comment=comment))

        fout.close()

        return legend_filename

    return None


# Templates for XML file (modified after Bruno Vellutini):
# Template for XML header.
BEGIN_TEMPLATE_1 = '''<?xml version="1.0" encoding="UTF-8"?>
<TrackMate version="3.4.2">
  <Model spatialunits="pixels" timeunits="frames">
    <FeatureDeclarations>
      <SpotFeatures>
        <Feature feature="QUALITY" name="Quality" shortname="Quality" dimension="QUALITY" isint="false" />
        <Feature feature="POSITION_X" name="X" shortname="X" dimension="POSITION" isint="false" />
        <Feature feature="POSITION_Y" name="Y" shortname="Y" dimension="POSITION" isint="false" />
        <Feature feature="POSITION_Z" name="Z" shortname="Z" dimension="POSITION" isint="false" />
        <Feature feature="POSITION_T" name="T" shortname="T" dimension="TIME" isint="false" />
        <Feature feature="FRAME" name="Frame" shortname="Frame" dimension="NONE" isint="true" />
        <Feature feature="RADIUS" name="Radius" shortname="R" dimension="LENGTH" isint="false" />
        <Feature feature="VISIBILITY" name="Visibility" shortname="Visibility" dimension="NONE" isint="true" />
        <Feature feature="SOURCE_ID" name="Source ID" shortname="Source" dimension="NONE" isint="true" />
        <Feature feature="CELL_DIVISION_TIME" name="Cell division time" shortname="Cell div. time" dimension="TIME" isint="false" />
        <Feature feature="MANUAL_COLOR" name="Manual spot color" shortname="Spot color" dimension="NONE" isint="true" />\n'''

BEGIN_METADATA_TEMPLATE = '''        <Feature feature="SHAPE" name="SIMIBioCell spot shape" shortname="Spot shape" dimension="NONE" isint="true" />
        <Feature feature="FATE" name="SIMIBioCell cell fate" shortname="Cell fate" dimension="NONE" isint="true" />
        <Feature feature="COMMENT" name="SIMIBioCell user comment" shortname="Comment" dimension="NONE" isint="true" />\n'''

BEGIN_TEMPLATE_2 = '''    </SpotFeatures>
      <EdgeFeatures>
        <Feature feature="SPOT_SOURCE_ID" name="Source spot ID" shortname="Source ID" dimension="NONE" isint="true" />
        <Feature feature="SPOT_TARGET_ID" name="Target spot ID" shortname="Target ID" dimension="NONE" isint="true" />
        <Feature feature="LINK_COST" name="Link cost" shortname="Cost" dimension="NONE" isint="false" />
        <Feature feature="VELOCITY" name="Velocity" shortname="V" dimension="VELOCITY" isint="false" />
        <Feature feature="DISPLACEMENT" name="Displacement" shortname="D" dimension="LENGTH" isint="false" />
        <Feature feature="MANUAL_COLOR" name="Manual edge color" shortname="Edge color" dimension="NONE" isint="true" />
      </EdgeFeatures>
      <TrackFeatures>
        <Feature feature="TRACK_INDEX" name="Track index" shortname="Index" dimension="NONE" isint="true" />
        <Feature feature="TRACK_ID" name="Track ID" shortname="ID" dimension="NONE" isint="true" />
        <Feature feature="TRACK_DURATION" name="Duration of track" shortname="Duration" dimension="TIME" isint="false" />
        <Feature feature="TRACK_START" name="Track start" shortname="T start" dimension="TIME" isint="false" />
        <Feature feature="TRACK_STOP" name="Track stop" shortname="T stop" dimension="TIME" isint="false" />
        <Feature feature="TRACK_DISPLACEMENT" name="Track displacement" shortname="Displacement" dimension="LENGTH" isint="false" />
        <Feature feature="NUMBER_SPOTS" name="Number of spots in track" shortname="N spots" dimension="NONE" isint="true" />
        <Feature feature="NUMBER_GAPS" name="Number of gaps" shortname="Gaps" dimension="NONE" isint="true" />
        <Feature feature="LONGEST_GAP" name="Longest gap" shortname="Longest gap" dimension="NONE" isint="true" />
        <Feature feature="NUMBER_SPLITS" name="Number of split events" shortname="Splits" dimension="NONE" isint="true" />
        <Feature feature="NUMBER_MERGES" name="Number of merge events" shortname="Merges" dimension="NONE" isint="true" />
        <Feature feature="NUMBER_COMPLEX" name="Complex points" shortname="Complex" dimension="NONE" isint="true" />
        <Feature feature="DIVISION_TIME_MEAN" name="Mean cell division time" shortname="Mean div. time" dimension="TIME" isint="false" />
        <Feature feature="DIVISION_TIME_STD" name="Std cell division time" shortname="Std div. time" dimension="TIME" isint="false" />
      </TrackFeatures>
    </FeatureDeclarations>\n'''

# Templates for spots section.
ALLSPOTS_BEGIN_TEMPLATE = '    <AllSpots nspots="{nspots}">\n'
INFRAME_BEGIN_TEMPLATE =  '      <SpotsInFrame frame="{frame}">\n'
SPOT_TEMPLATE =           '        <Spot ID="{id}" name="{name}" VISIBILITY="1" RADIUS="{radius}" QUALITY="-1.0" SOURCE_ID="0" POSITION_T="{t}" POSITION_X="{x}" POSITION_Y="{y}" FRAME="{frame}" POSITION_Z="{z}" '
SPOT_COLOR_TEMPLATE =     'MANUAL_COLOR="{color}" '
SPOT_SHAPE_TEMPLATE =     'SHAPE="{shape}" '
SPOT_FATE_TEMPLATE =      'FATE="{fate}" '
SPOT_COMMENT_TEMPLATE =   'COMMENT="{id}" '
INFRAME_END_TEMPLATE =    '      </SpotsInFrame>\n'
ALLSPOTS_END_TEMPLATE =   '    </AllSpots>\n'

# Templates for tracks section.
ALLTRACKS_BEGIN_TEMPLATE = '    <AllTracks>\n'
TRACK_BEGIN_TEMPLATE =     '      <Track name="{name}" TRACK_INDEX="{id}" TRACK_ID="{id}" TRACK_DURATION="{duration}" TRACK_START="{start}" TRACK_STOP="{stop}" TRACK_DISPLACEMENT="{displacement}" NUMBER_SPOTS="{nspots}" NUMBER_GAPS="{ngaps}" LONGEST_GAP="{longestgap}" NUMBER_SPLITS="{nsplits}" NUMBER_MERGES="0" NUMBER_COMPLEX="0" DIVISION_TIME_MEAN="NaN" DIVISION_TIME_STD="NaN">\n'
EDGE_TEMPLATE =            '        <Edge SPOT_SOURCE_ID="{source}" SPOT_TARGET_ID="{target}" LINK_COST="-1.0" VELOCITY="{velocity}" DISPLACEMENT="{displacement}" '
EDGE_COLOR_TEMPLATE =      'MANUAL_COLOR="{color}" '
TRACK_END_TEMPLATE =       '      </Track>\n'
ALLTRACKS_END_TEMPLATE =   '    </AllTracks>\n'

# Templates for filtered tracks section.
FILTEREDTRACKS_BEGIN_TEMPLATE = '    <FilteredTracks>\n'
FILTEREDTRACKS_TEMPLATE =       '      <TrackID TRACK_ID="{id}" />\n'
FILTEREDTRACKS_END_TEMPLATE =   '    </FilteredTracks>\n'

# Templates for SIMI BioCell legends section.
BIOCELLLEGENDS_BEGIN_TEMPLATE = '    <SIMIBioCellLegends>\n'
SPOT_SHAPES_TEMPLATE =          '      <SpotShapes Sphere="0" />\n      <SpotShapes Cube="1" />\n      <SpotShapes Spiky="2" />\n      <SpotShapes Diamond="3" />\n      <SpotShapes Ufo="4" />\n'
CELL_FATES_TEMPLATE =           '      <CellFates {fate}="{index}" />\n'
USER_COMMENT_TEMPLATE =         '      <UserComment Spot_ID{id}="{comment}" />\n'
BIOCELLLEGENDS_END_TEMPLATE =   '    </SIMIBioCellLegends>\n'

# Template for image data.
IMAGE_DATA_TEMPLATE = '  </Model>\n  <Settings>\n    <ImageData filename="{file}" folder="{dir}" width="{width}" height="{height}" nslices="{nslices}" nframes="{nframes}" pixelwidth="{pixelwidth}" pixelheight="{pixelheight}" voxeldepth="{voxeldepth}" timeinterval="{timeinterval}" />\n'

# Template for XML ending.
END_TEMPLATE = '''    <InitialSpotFilter feature="QUALITY" value="0.0" isabove="true" />
    <SpotFilterCollection />
    <TrackFilterCollection />
    <AnalyzerCollection>
      <SpotAnalyzers>
        <Analyzer key="Spot Source ID" />
        <Analyzer key="CELL_DIVISION_TIME_ON_SPOTS" />
        <Analyzer key="MANUAL_SPOT_COLOR_ANALYZER" />
      </SpotAnalyzers>
      <EdgeAnalyzers>
        <Analyzer key="Edge target" />
        <Analyzer key="Edge velocity" />
        <Analyzer key="MANUAL_EDGE_COLOR_ANALYZER" />
      </EdgeAnalyzers>
      <TrackAnalyzers>
        <Analyzer key="Track index" />
        <Analyzer key="Track duration" />
        <Analyzer key="Branching analyzer" />
        <Analyzer key="CELL_DIVISION_TIME_ANALYZER" />
      </TrackAnalyzers>
    </AnalyzerCollection>
  </Settings>
  <GUIState>
    <SetupAssignments>
      <ConverterSetups>
        <ConverterSetup>
          <id>0</id>
          <min>90.0</min>
          <max>220.0</max>
          <color>-1</color>
          <groupId>0</groupId>
        </ConverterSetup>
      </ConverterSetups>
      <MinMaxGroups>
        <MinMaxGroup>
          <id>0</id>
          <fullRangeMin>0.0</fullRangeMin>
          <fullRangeMax>65535.0</fullRangeMax>
          <rangeMin>0.0</rangeMin>
          <rangeMax>65535.0</rangeMax>
          <currentMin>90.0</currentMin>
          <currentMax>220.0</currentMax>
        </MinMaxGroup>
      </MinMaxGroups>
    </SetupAssignments>
    <Bookmarks />
  </GUIState>
</TrackMate>'''

# Templates for legend.txt file:
TITLE_TEMPLATE =          'LEGENDS FOR FILE: {file}\n\n'
ODD_MITOSES_HEAD =        'ODD MITOSES (Polytomies and mitoses with only one offspring cell) - Consider revising:\n\n'
ONE_OFFSPRING_HEAD =      'Mitoses with only one offspring cell:\n'
MULTIPLE_OFFSPRING_HEAD = 'Polytomies (i.e. mitoses include alternative branch(es) in SIMI BioCell lineage):\n'
ODD_MITOSES_TEMPLATE =    '  Parent Cell: {name}\n'
SHAPES_TEMPLATE =         'SPOT SHAPES:\nIndex\tShape\n0\t-\tSphere\n1\t-\tCube\n2\t-\tSpiky\n3\t-\tDiamond\n4\t-\tUfo\n\n'
FATES_HEAD_TEMPLATE =     'CELL FATES:\nIndex\tFate\n'
FATES_TEMPLATE =          '{index}\t-\t{fate}\n'
COMMENT_HEAD_TEMPLATE =   'USER COMMENTS:\nSpotID\t\t\tComment\n'
COMMENT_TEMPLATE =        '{id}\t\t-\t\t{comment}\n'

if __name__ == '__main__':
    biocell2xml()
