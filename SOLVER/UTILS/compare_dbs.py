#!/usr/bin/env python
# coding: utf-8
import argparse
import instaseis
import matplotlib.pyplot as plt


def define_arguments():
    helptext = 'Compare two databases.'
    formatter_class = argparse.RawTextHelpFormatter
    parser = argparse.ArgumentParser(description=helptext,
                                     formatter_class=formatter_class)

    helptext = "Database 1 directory name. \n"
    parser.add_argument('db_1_name', help=helptext)
    helptext = "Database 2 directory name. \n"
    parser.add_argument('db_2_name', help=helptext)

    return parser


if __name__ == "__main__":
    parser = define_arguments()
    args = parser.parse_args()

    db1 = instaseis.open_db(args.db_1_name)
    db2 = instaseis.open_db(args.db_2_name)

    src = instaseis.Source(latitude=0.0, longitude=0.0, m_rr=1e20)

    rec = instaseis.Receiver(latitude=45.0, longitude=0.0)

    st1 = db1.get_seismograms(src, rec)
    st2 = db2.get_seismograms(src, rec)

    plt.plot(st1[0].times(), st1[0].data, label=args.db_1_name)
    plt.plot(st2[0].times(), st2[0].data, label=args.db_2_name)

    plt.xlabel('time / seconds')
    plt.legend()

    plt.show()
