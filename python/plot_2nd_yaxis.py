import sys

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt  # noqa


def plot_2nd_yaxis(x, y, z, fmt1, fmt2,
                   xlab, ylab, label1, label2, output):
    fig, x1 = plt.subplots()
    lines = list()
    line, = x1.plot(x, y, fmt1)
    lines.append(line)
    x1.set_xlabel(xlab)
    x1.set_ylabel(ylab)
    line, = x1.plot(x, z, fmt2)
    lines.append(line)
    x1.legend(lines, (label1, label2), loc='best')
    fig.tight_layout()
    plt.savefig(output)


if __name__ == "__main__":
    if len(sys.argv) < 8:
        sys.exit(('Usage: plot_2nd_yaxis.py infile output xlab ylab '
                  'label1 label2 ratio'))
    infile, output, xlab, ylab, label1, label2, ratio = sys.argv[1:]
    x = list()
    y = list()
    z = list()
    with open(infile) as f:
        f.readline()
        for row in f:
            rx, ry, rz = row.strip().split('\t')
            x.append(float(rx))
            y.append(float(ry))
            z.append(float(rz) * float(ratio))
    plot_2nd_yaxis(x, y, z, 'r-', 'b--', xlab, ylab, label1, label2, output)
