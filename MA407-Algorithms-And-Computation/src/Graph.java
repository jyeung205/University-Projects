// Candidate Number 36957
//--------------------------------------------------------------

public class Graph {
    Node nodeListHead;
    Node nodeListEnd;
    Edge edgeListHead;
    Edge edgeListEnd;
    final int maxCoordinate = 100;
    final double textWidth = 17.0;

    public static class Node {
        int x;
        int y;
        int label;
        String value;
        Node next;

        Node(int x, int y, int label, String value) {
            this.x = x;
            this.y = y;
            this.label = label;
            this.value = value;
        }
    }

    public static class Edge {
        int nodeLabel1;
        int nodeLabel2;
        Edge next;
    }

    void addNode(int x, int y, int label, String value) {

        // Check to see if coordinates exceed maxCoordinate value
        if (Math.abs(x) > maxCoordinate && Math.abs(y) > maxCoordinate) {
            System.out.println("% coordinate " + x + " and " + y + " in node (" + x + "," + y + ") is too large");
            return;
        }
        if (Math.abs(x) > maxCoordinate) {
            System.out.println("% coordinate " + x + " in node (" + x + "," + y + ") is too large");
            return;
        }
        if (Math.abs(y) > maxCoordinate) {
            System.out.println("% coordinate " + y + " in node (" + x + "," + y + ") is too large");
            return;
        }

        // If the nodeList is empty then add a new node to the list
        if (nodeListHead == null) {
            nodeListHead = new Node(x, y, label, value);
            nodeListEnd = nodeListHead;
        }
        else {
            // Check for repeated coordinates and update value if coordinates are repeated
            for (Node node = nodeListHead; node != null; node = node.next) {
                if (node.x == x && node.y == y) {
                    node.value = value;
                    return;
                }
            }
            // Else create a new node and add to the list
            Node node = new Node(x, y, label, value);
            nodeListEnd.next = node;
            nodeListEnd = node;
        }
    }

    void addEdge(int x1, int y1, int x2, int y2) {

        // find label of x1, y1 and x2, y2
        Edge edge = new Edge();
        boolean node1Present = false;
        boolean node2Present = false;
        for (Node node = nodeListHead; node != null; node = node.next) {
            if (node.x == x1 && node.y == y1) {
                edge.nodeLabel1 = node.label;
                node1Present = true;
            }
            if (node.x == x2 && node.y == y2) {
                edge.nodeLabel2 = node.label;
                node2Present = true;
            }
        }

        // Check to see if nodes connected by the edge are present before adding them to the edgeList
        if (!node1Present && !node2Present) {
            System.out.println("% Tried to add edge from (" + x1 + "," + y1 + ") to (" + x2 + "," + y2 + ") but both nodes are not present in the graph");
            return;
        }
        if (!node1Present) {
            System.out.println("% Tried to add edge from (" + x1 + "," + y1 + ") to (" + x2 + "," + y2 + ") but node with coordinates (" + x1 + "," + y1 + ") is not present in the graph");
            return;
        }
        if (!node2Present) {
            System.out.println("% Tried to add edge from (" + x1 + "," + y1 + ") to (" + x2 + "," + y2 + ") but node with coordinates (" + x2 + "," + y2 + ") is not present in the graph");
            return;
        }

        if (edgeListHead == null) {
            edgeListHead = edge;
        }
        else {
            edgeListEnd.next = edge;
        }
        edgeListEnd = edge;

    }

    void clear() {
        nodeListHead = null;
        nodeListEnd = null;
        edgeListHead = null;
        edgeListEnd = null;
    }

    void outHeader() {
        double oddSideMargin = -1* (2.54 - (21.0 - textWidth) / 2);

        System.out.print("\\documentclass[a4paper,11pt]{article}\n" +
                "\\usepackage{mathpazo}\n" +
                "\\usepackage{tikz}\n" +
                "\\usetikzlibrary{shapes}\n");

        System.out.print("\\oddsidemargin "); System.out.printf("%.2f", oddSideMargin); System.out.print("cm\n");

        System.out.print("\\textwidth "); System.out.printf("%.2f", textWidth); System.out.print("cm\n");

        System.out.println("\\textheight 24cm\n" +
                "\\topmargin -1.3cm\n" +
                "\\parindent 0pt\n" +
                "\\parskip 1ex\n" +
                "\\pagestyle{empty}\n" +
                "\\begin{document}\n" +
                "\\medskip\\hrule\\medskip\n");

    }

    void outFooter() {
        System.out.println("\\medskip\\hrule\\medskip\n" +
                "\\end{document}");
    }

    void outGraph() {

        // Loop through nodeList and output LaTex command
        for (Node node = nodeListHead; node != null; node = node.next) {
            int x = node.x;
            int y = node.y;
            int label = node.label;
            String value = node.value;

            System.out.print( "\\draw [thick] "+ "(");
            System.out.print(x + "," + y);
            System.out.print(") node[draw, rounded rectangle] ");
            System.out.println("(" + label + ")" + " {" + value + "};");
        }

        // Loop through edgeList and output LaTex command
        for (Edge edge = edgeListHead; edge != null; edge = edge.next) {
            int label1 = edge.nodeLabel1;
            int label2 = edge.nodeLabel2;
            System.out.println("\\draw [->, thick] (" + label1 + ") to (" + label2 + ");");
        }

    }

    // Find the smallest and largest x, y coordinates and return an array containing them
    int[] getGridCoordinates() {

        int xMax = -100;
        int yMax = -100;
        int xMin = 100;
        int yMin = 100;
        for (Node node = nodeListHead; node != null; node = node.next) {
            if (node.x < xMin && node.x >= -maxCoordinate) {
                xMin = node.x;
            }
            if (node.x > xMax && node.x <= maxCoordinate) {
                xMax = node.x;
            }
            if (node.y < yMin && node.y >= -maxCoordinate) {
                yMin = node.y;
            }
            if (node.y > yMax && node.y <= maxCoordinate) {
                yMax = node.y;
            }
        }

        int[] gridCoordinates = new int[4];
        gridCoordinates[0] = xMin;
        gridCoordinates[1] = yMin;
        gridCoordinates[2] = xMax;
        gridCoordinates[3] = yMax;
        return gridCoordinates;
    }

    // Outputs the LaTex commands to draw a grid between the smallest and largest coordinates
    void outGrid(int[] gridCoordinates) {

        int xMin = gridCoordinates[0];
        int yMin = gridCoordinates[1];
        int xMax = gridCoordinates[2];
        int yMax = gridCoordinates[3];
        System.out.println("\\draw [help lines, color=green] (" + xMin + "," + yMin + ") grid " + "(" + xMax + "," + yMax + ");\n");

    }

    // Method to output LaTex commands
    void latex(String[] args, int[] gridCoordinates, Tree tree) {

        if (args != null) {
            for (String arg : args) {
                System.out.print(arg + " ");
            }
            System.out.println("are inserted \n");
        }

        // printSorted and averageDepth are called when drawing a tree but not a graph
        if (tree != null) {
            System.out.print("In sorted order: ");
            tree.printSorted(tree.root);
            System.out.println("\n");

            tree.createDepthList(tree.root, 0);
            tree.averageDepth(tree.depthListHead, true);
        }

        double scale = 0.600;
        int xMin = gridCoordinates[0];
        int xMax = gridCoordinates[2];
        int width = 1;

        // Calculate width of grid
        if (xMin > 0 && xMax > 0 || xMax > 0 && xMin < 0) {
            width = xMax - xMin + 1;
        }
        if (xMin < 0 && xMax < 0) {
            width = Math.abs(xMin) - Math.abs(xMax) + 1;
        }

        // 28 is the max width for textWidth 17.0cm and scale 0.600cm
        if (width > 28) {
            scale = textWidth / width;
        }
        if (scale < 0.3) {
            scale = 0.3;
        }

        // Output LaTex commands to draw grid and graph
        System.out.print("\\begin{tikzpicture}");
        System.out.printf("[scale=%.3f]", scale);
        System.out.println();
        outGrid(gridCoordinates);
        outGraph();
        System.out.println("\n\\end{tikzpicture}\n");
    }

    public static void main(String[] args) {
	    Graph graph = new Graph();

	    // Take the command line inputs and add the nodes in the graph
        for (int i=0; i<args.length/2; i++) {
            int x = Integer.parseInt(args[2*i]);
            int y = Integer.parseInt(args[2*i+1]);
            graph.addNode(x, y, i, Integer.toString(i+1));
        }

        // Take the command line inputs and add the edges in the graph
        for (int i=0; i<args.length/2 - 1; i++) {
            int x1 = Integer.parseInt(args[2*i]);
            int y1 = Integer.parseInt(args[2*i+1]);
            int x2 = Integer.parseInt(args[2*i+2]);
            int y2 = Integer.parseInt(args[2*i+3]);
            graph.addEdge(x1, y1, x2, y2);
        }

        // Call methods to output LaTex commands
        graph.outHeader();

        // If no input arguments
        if (args.length == 0) {
            System.out.println("Please input some coordinates \n");
        }

        int[] gridCoordinates = new int[4];

        // If no nodes were added to the graph
        if (graph.nodeListHead == null) {
            System.out.println("The graph is empty");
        }
        else {
            gridCoordinates = graph.getGridCoordinates();
        }

        graph.latex(null, gridCoordinates, null);
        graph.outFooter();

    }
}
