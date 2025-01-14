#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt

show_id = 1
show_id_only_pin = 1

parser = argparse.ArgumentParser(description='Draw Nets/Trees')
parser.add_argument('--input_filenames', type=str, default=[
    "./toys/DelayTree_toy1.tree",
    ])

parser.add_argument('-a', '--anno_text', choices=['no', 'pin_only', 'all'], default='all')
parser.add_argument('-s', '--scale', type=float, default=1.0)
args = parser.parse_args()

types = ['Net', 'Tree']

class Node:
    def __init__(self, line, withParent = False):
        tokens = line.split()
        self.id, self.x, self.y = int(tokens[0]), float(tokens[1]), float(tokens[2])
        if withParent:
            self.parent_id = int(tokens[3])

for filename in args.input_filenames:
    with open(filename) as file:
        type = None
        line = file.readline()
        while line != '':
            for t in types:
                if t in line:
                    type = t
                    break
            if type is not None:
                break
            line = file.readline()
        
        if type == 'Net' or type == 'Tree':
            # header
            tokens = line.split()
            net_id, net_name, num_pins = int(tokens[1]), tokens[2], int(tokens[3])
            nodes = []
            line = file.readline()
            while line != '':
                nodes.append(Node(line, type == 'Tree')) # (id, x, y, parent_id)
                line = file.readline()
            plt.figure(figsize=(15, 15))

            # edges
            if type == 'Tree':
                for node in nodes:
                    if node.id == 0:
                        continue
                    parent = nodes[node.parent_id]
                    plt.plot([node.x, parent.x], [node.y, parent.y], c='k')

            labeled_node = [1,2,10,15,18,20]
            # nodes
            for node in nodes:
                if node.id == 0: # source
                    plt.plot(node.x, node.y, c='r', marker='s', mew=0, ms=6.5 * args.scale)
                elif node.id < num_pins: # sinks
                    # plt.plot(node.x, node.y, c='k', marker='s', mew=0, ms=6.5 * args.scale)
                    if node.id in labeled_node:
                        plt.plot(node.x, node.y, c='r', marker='*', mew=0, ms=15.5 * args.scale)
                    else:
                        plt.plot(node.x, node.y, c='k', marker='s', mew=0, ms=6.5 * args.scale)
                # else: # Steiner points
                #     plt.plot(node.x, node.y, c='b', marker='o', mew=0, ms=15.5 * args.scale)
                if show_id :
                    if show_id_only_pin:
                        if node.id < num_pins:
                            if node.id in labeled_node:
                                plt.text(node.x+0.66, node.y+0.66, node.id, fontsize=40, color='red',)
                            else:
                                plt.text(node.x+0.66, node.y+0.66, node.id, fontsize=40)
                    else:
                        plt.text(node.x+0.66, node.y+0.66, node.id, fontsize=40)
                    # plt.text(node.x, node.y, node.id, fontsize=6)
                # anno text
                if args.anno_text == 'no':
                    continue
                if node.id < 0:
                    continue
                if args.anno_text == 'pin_only' and node.id >= num_pins:
                    continue
                # if node.id < num_pins: # source
                #     plt.text(node.x + .5, node.y, str(node.id), c='r', fontsize=7)
                # else:
                #     plt.text(node.x - .5, node.y + .5, str(node.id), c='b', fontsize=7)

            # format & save
            plt.axis('square')
            plt.axis('off')
            plt.savefig('{}.png'.format(filename), bbox_inches='tight')
            # Moving files to another folder
            import shutil
            move_filename = filename.split('/')[-1]
            shutil.move('{}.png'.format(filename), './run/{}.png'.format(move_filename))

        else:
            print("Unknown type: ", type)