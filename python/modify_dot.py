#!/bin/env python
import re, os
from collections import defaultdict
import imageio
import argparse, textwrap

class dotgraph:
    def __init__(self, dotfile):
        self.dir, self.dotfile = os.path.split(dotfile)
        self.lines = []
        self.colormap = {}
        self.read_dot(dotfile)

    def read_dot(self, filename):
        self.lines = []
        with open(filename) as f:
            for line in f:
                while line.count("[") != line.count("]"):
                    line = line+f.next()
                self.lines.append(line)
                m = re.search("original_id=\"\[?u?\'?(.*?)(\'\])?\"", line)
                if m:
                    id = m.group(1)
                    m = re.search("fillcolor=\"(#\w*)\"", line)
                    if m:
                        color = m.group(1)
                        self.update_colormap(id, color)


    def write_dot(self, outfile):
        with open(outfile, 'w') as o:
            for line in self.lines:
                o.write(line)

    def change_color(self):
        lines = []
        for line in self.lines:
            m = re.search("original_id=\"\[?u?\'?(.*?)(\'\])?\"", line)
            if m:
                id = m.group(1)
                if id in self.colormap:
                    line = re.sub("fillcolor=\"#\w*\"", "fillcolor=\"{}\"".format(self.colormap[id]), line)
            lines.append(line)
        self.lines = lines

    def update_colormap(self, id, color):
        self.colormap[id] = color

    def read_colormap(self, mapfile, ididx=0, coloridx=1):
        with open(mapfile) as f:
            for line in f:
                elements = line.strip().split()
                self.update_colormap(elements[ididx], elements[coloridx])

def make_step_graph(dotgraph, output, newcolormap):
    colorlist = sorted(set(dotgraph.colormap.values()), reverse=True) #lighter color appears earlier
    color_to_id = defaultdict(list)
    graphfilenames = []
    os.mkdir(output)

    for id, color in dotgraph.colormap.iteritems():
        color_to_id[color].append(id)

    if newcolormap:
        dotgraph.read_colormap(newcolormap)

    for i, target in enumerate(colorlist):
        for id in color_to_id[target]: #change target color to white
            dotgraph.update_colormap(id, "#ffffff") #change color to white
        dotgraph.change_color()
        filename = dotgraph.dotfile+"."+"%03d" % (len(colorlist)-i)+".dot" #let the number starts from the largest one
        dotgraph.write_dot(os.path.join(output, filename))
        graphfilenames.append(os.path.join(output, filename+".png"))

    # os.system("dot -Tpng -s300 movie/*.dot -O")
    #
    # with imageio.get_writer(dotgraph.dotfile+".gif", mode='I', duration=0.5) as writer:
    #     for filename in graphfilenames:
    #         image = imageio.imread(filename)
    #         writer.append_data(image)

if __name__ == '__main__':
    #read the arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Read dot file then change node color based on colormap file.
        Multiple stepwise dot files could be generated automatically with --make-step-graph set to True.
        '''))
    parser.add_argument('-s', '--make-step-graph', action='store_true', help='Generate stepwise dot files into output path')
    parser.add_argument('-i', '--input', type=str, help='Input file name', required=True)
    parser.add_argument('-m', '--colormap', type=str, help='Colormap files used to change the node color')
    parser.add_argument('-o', '--output', type=str, help='Output file name or path', required=True)
    args = parser.parse_args()

    gp = dotgraph(args.input)

    if args.make_step_graph:
        make_step_graph(gp, args.output, args.colormap)
    else:
        if args.colormap:
            gp.read_colormap(args.colormap)
            gp.change_color()
        gp.write_dot(args.output)
