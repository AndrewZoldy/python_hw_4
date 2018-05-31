from Bio import SeqIO
from graphviz import Digraph
import argparse
import os

class Vertex:

    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.inner_edges = {}
        self.outer_edges = {}

    def increase_coverage(self):
        self.coverage += 1


class Edge:

    def __init__(self, k1, k2):
        self.seq = k1 + k2[-1]
        self.coverage = 0
        self.number = 2

    def calc_coverage(self, c1, c2):
        self.coverage = c1 + c2


class Graph:

    def __init__(self, size):
        self.vertices = {}
        self.s = size
        self.graph = Digraph()

    def add_read(self, read):
        read_lng = len(read)
        if read_lng < self.s:
            return

        kmer = read[:self.s]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)

        for next_kmer_indx in range(1, read_lng - self.s + 1, 1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx + self.s)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)

            new_edge = Edge(kmer, next_kmer[-1])

            self.vertices[next_kmer].inner_edges[kmer] = new_edge

            self.vertices[kmer].outer_edges[next_kmer] = new_edge

            kmer = next_kmer

    def calc_init_edge_coverage(self):

        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].outer_edges.keys():
                self.vertices[current_vertex].outer_edges[next_vertex].calc_coverage(
                    self.vertices[current_vertex].coverage, self.vertices[next_vertex].coverage)

    def graphviz(self, output, mode):

        for vertex in self.vertices:
            if mode == 'full':
                self.graph.node(vertex, label=vertex)
                for next_vertex in self.vertices[vertex].outer_edges:
                    self.graph.edge(vertex, next_vertex, label=self.vertices[vertex].outer_edges[next_vertex].seq)
            elif mode == 'cutted':
                self.graph.node(vertex, label=str(self.vertices[vertex].coverage))
                for next_vertex in self.vertices[vertex].outer_edges:
                    self.graph.edge(vertex, next_vertex,
                                    label=str(self.vertices[vertex].outer_edges[next_vertex].coverage) + ',' +
                                          str(len(self.vertices[vertex].outer_edges[next_vertex].seq)))

        self.graph.render(filename=output, view=True)

    def collapse_graph(self):

        needed_vertices = []
        for vertex in self.vertices:
            if len(self.vertices[vertex].inner_edges) == 1 and len(self.vertices[vertex].outer_edges) == 1:
                needed_vertices.append(self.vertices[vertex].seq)

        for i in range(len(needed_vertices)):
            if needed_vertices[i] in self.vertices and len(self.vertices) > 2:
                inner_vertex = [j for j in self.vertices[needed_vertices[i]].inner_edges.keys()][0]
                outer_vertex = [m for m in self.vertices[needed_vertices[i]].outer_edges.keys()][0]
                new_edge = Edge(self.vertices[needed_vertices[i]].inner_edges[inner_vertex].seq,
                                self.vertices[needed_vertices[i]].outer_edges[outer_vertex].seq[-1])
                self.vertices[inner_vertex].outer_edges[outer_vertex] = new_edge
                self.vertices[outer_vertex].inner_edges[inner_vertex] = new_edge
                # edge length
                self.vertices[inner_vertex].outer_edges[outer_vertex].n = self.vertices[outer_vertex].inner_edges[
                                                                       needed_vertices[i]].n + \
                                                                   self.vertices[inner_vertex].outer_edges[
                                                                       needed_vertices[i]].n - 1
                self.vertices[outer_vertex].inner_edges[inner_vertex].n = self.vertices[outer_vertex].inner_edges[
                                                                      needed_vertices[i]].n + \
                                                                  self.vertices[inner_vertex].outer_edges[
                                                                      needed_vertices[i]].n - 1
                # edge coverage
                self.vertices[inner_vertex].outer_edges[outer_vertex].coverage = self.vertices[inner_vertex].outer_edges[
                                                                              needed_vertices[i]].coverage + \
                                                                          self.vertices[outer_vertex].coverage
                self.vertices[outer_vertex].inner_edges[inner_vertex].coverage = self.vertices[inner_vertex].outer_edges[
                                                                             needed_vertices[i]].coverage + \
                                                                         self.vertices[outer_vertex].coverage

                del self.vertices[outer_vertex].inner_edges[needed_vertices[i]]
                del self.vertices[inner_vertex].outer_edges[needed_vertices[i]]
                del self.vertices[needed_vertices[i]]
            else:
                continue

    def obtained_edge_coverage(self):
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].outer_edges.keys():
                self.vertices[current_vertex].outer_edges[next_vertex].coverage = self.vertices[current_vertex].outer_edges[
                                                                                    next_vertex].coverage / \
                                                                                self.vertices[current_vertex].outer_edges[
                                                                                    next_vertex].n

    def write_out(self, fasta, size):
        id = 0
        with open(fasta, 'w') as f:
            for vertex in self.vertices:
                for next_vertex in self.vertices[vertex].outer_edges:
                    if len(self.vertices[vertex].outer_edges[next_vertex].seq) > size + 1:
                        id += 1
                        f.write('>'+'id_'+str(id)+'\n')
                        f.write(str(self.vertices[vertex].outer_edges[next_vertex].seq)+'\n')
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Graph')

    parser.add_argument('-i', '--input', help='path to input fasta file', type=str, required=True)
    parser.add_argument('-o', '--output_graph', help='path to output .dot file with graph', type=str, required=True)
    parser.add_argument('-s', '--size', help='kmer size', type=int, default=15)
    parser.add_argument('-m', '--mode', help='choose full or cutted view', type=str, default='full')
    parser.add_argument('-c', '--collapse', help='switch to False if you don not want to collapse graph', action='store_true', default=True)
    parser.add_argument('-of', '--output_fasta', help='path to output .fasta file wich will include assembly', type=str, default=os.curdir+'/output_graph_fasta.fasta')

    args = parser.parse_args()
    file = args.input
    collapse = args.collapse
    size = args.size
    mode = args.mode
    output = args.output_graph
    fasta=args.output_fasta


    graph_exemplar = Graph(size)

    with open(file, 'r') as input:
        for record in SeqIO.parse(input, 'fasta'):
            read = str(record.seq)
            graph_exemplar.add_read(read)
            graph_exemplar.add_read(str(record.reverse_complement().seq))


    if collapse == True:
        graph_exemplar.collapse_graph()
    else:
        pass

    graph_exemplar.obtained_edge_coverage()
    graph_exemplar.graphviz(output, mode)

    if collapse == True:
        graph_exemplar.write_out(fasta, size)