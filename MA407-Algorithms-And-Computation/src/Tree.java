// Candidate Number 36957
//-----------------------------------------------------

public class Tree {
    Item root = null;
    DepthItem depthListHead;
    DepthItem depthListEnd;

    public static class Item {  // Item for binary search tree
        int  value, x, y;
        Item left, right, parent;

        private Item (int val, Item p) {
            value = val;
            parent = p;
        }
    }

    public static class DepthItem { // Item for calculating average depth of tree
        int depth;
        int count;
        DepthItem next;

        DepthItem(int depth, int count) {
            this.depth = depth;
            this.count = count;
        }
    }

    // find method given in BinarySearchTree.java file
    Item find(int x) {
        Item p = root;
        while (p != null && p.value != x) {
            if (x < p.value) {
                p = p.left;
            }
            else {
                p = p.right;
            }
        }
        return p;
    }

    // insert method given in BinarySearchTree.java file
    Item insert(int x) {
        if (root == null) {
            root = new Item(x, null);
            return root;
        }
        Item p = root, q = null;
        while (p != null) {
            q = p;
            if (x == p.value) { // value found
                return p;
            }
            else if (x < p.value) {
                p = p.left;
            }
            else {
                p = p.right;
            }
        }
        p = new Item(x, q);
        if (x < q.value) {
            q.left = p;
        }
        else {
            q.right = p;
        }
        return p;
    }

    // delete method given in BinarySearchTree.java file
    void delete(int x) {
        Item p = find(x);
        if (p == null) {
            return;
        }
        else if (p.right == null) {
            deleteWithoutRight(p);
        }
        else if (p.left == null) {
            deleteWithoutLeft(p);
        }
        else {
            Item q = p.right;
            while (q.left != null)
                q = q.left;
            p.value = q.value;
            deleteWithoutLeft(q);
        }
    }

    // deleteWithoutLeft method given in BinarySearchTree.java file
    void deleteWithoutLeft(Item p)  {
        Item q = p.parent;
        if (q == null) {
            root = p.right;
            if (root != null) {
                root.parent = null;
            }
        }
        else if (q.left == p) {
            q.left = p.right;
            if (q.left != null) {
                q.left.parent = q;
            }
        }
        else {
            q.right = p.right;
            if (q.right != null) {
                q.right.parent = q;
            }
        }
    }

    // deleteWithoutRight method given in BinarySearchTree.java file
    void deleteWithoutRight(Item p) {
        Item q = p.parent;
        if (q == null) {
            root = p.left;
            if (root != null) {
                root.parent = null;
            }
        }
        else if (q.left == p) {
            q.left = p.left;
            if (q.left != null) {
                q.left.parent = q;
            }
        }
        else {
            q.right = p.left;
            if (q.right != null) {
                q.right.parent = q;
            }
        }
    }


    void printSorted(Item root) {
        if (root != null) {
            printSorted(root.left);
            System.out.print(root.value + " ");
            printSorted(root.right);
        }
    }

    // Method to create depthList which is used to calculate average depth
    void createDepthList(Item root, int depth) {

        if (depthListHead == null) {
            depthListHead = new DepthItem(depth, 1);
            depthListEnd = depthListHead;
        }

        if (root != null) {
            createDepthList(root.left, depth + 1);
            if (depth != 0) {

                // Check if depth level is already in the list. If it is then add one to the count
                boolean present = false;
                for (DepthItem depthItem = depthListHead; depthItem != null; depthItem = depthItem.next) {
                    if (depthItem.depth == depth) {
                        depthItem.count++;
                        present = true;
                    }
                }

                // Else add a new depthItem with count 1
                if (!present) {
                    DepthItem depthItem = new DepthItem(depth, 1);
                    depthListEnd.next = depthItem;
                    depthListEnd = depthItem;
                }
            }
            createDepthList(root.right, depth + 1);
        }
    }

    // Calculates the average depth of a tree
    double averageDepth(DepthItem depthListHead, boolean print) {
        int totalCount = 0;
        int sum = 0;
        for (DepthItem depthItem = depthListHead; depthItem != null; depthItem = depthItem.next) {
            sum += depthItem.depth * depthItem.count;
            totalCount += depthItem.count;
        }
        double averageDepth = (double) sum / totalCount;

        if (print){
            System.out.print("Average depth = ");
            System.out.printf("%.3f", averageDepth);
            System.out.println(", size " + totalCount + "\n");
        }

        return averageDepth;
    }

    // Creates a random permutation of numbers 1 to n
    String[] randomPerm(int n) {
        // -------------------------------------------
        // Argument for correctness:
        // We need to prove that the probability that the j th element goes to the k th position is 1/n for each position k in the array and j arbitrary.
        //
        // When k=n-1 (last position), this is clear as k is swapped with an index randomly chosen between 0 and n-1
        //
        // When k=n-2, there are two cases.
        //      Case 1: j = n-1 (last element)
        //      (Probability that the last element doesn't stay at the last position in the previous iteration) x (Probability the new index of the last element is selected in the current iteration) = (n-1) / n * 1 / (n-1) = 1/n
        //
        //      Case 2: 0<= j < n-1
        //      (Probability index j is not picked in previous iteration) x (Probability index j is picked in current iteration) = (n-1)/n * 1/(n-1) = 1/n
        //
        // Continuing with the above argument, we can show the probability the j th element ends up in position k, for k=n-1, n-2, ... , 0 is 1/n.
        // -------------------------------------------

        // Create an array with numbers 1 to n in ascending order
        String[] arr = new String[n];
        for (int i = 0; i < n; i++) {
            arr[i] = Integer.toString(i+1);
        }

        // Fisher–Yates shuffle
        for (int i=n-1; i>0; i--) {
            int random = (int) (Math.random() * (i+1));
            String swap = arr[random];
            arr[random] = arr[i];
            arr[i] = swap;
        }

        return arr;
    }

    // Generates 20 random permutations for n=3
    void randomTest() {
        for (int i=0; i < 20; i++) {
            String[] arr = randomPerm(3);
            for (String a : arr) {
                System.out.print(a + " ");
            }
            System.out.println();
        }
    }

    // Find and assign coordinates of the nodes in the tree
    int x = -100; // Left most x coordinate is set to -100
    void findCoordinates(Item item, int y) {
        if (item != null) {
            findCoordinates(item.left, y-2);
            item.x = x;
            item.y = y;
            x++;
            findCoordinates(item.right, y-2);
        }
    }

    // Adds nodes to graph class
    int label = 0; // Label for nodes in Graph class
    void treeGraphAddNodes(Item item, Graph graph) {
        if (item != null) {
            treeGraphAddNodes(item.left, graph);
            graph.addNode(item.x, item.y, label, Integer.toString(item.value));
            label++;
            treeGraphAddNodes(item.right, graph);
        }
    }

    // Add edges to graph class
    void treeGraphAddEdges(Item item, Graph graph) {
        if (item != null) {
            if (item.left != null) {
                treeGraphAddEdges(item.left, graph);
                graph.addEdge(item.x, item.y, item.left.x, item.left.y);
            }
            if (item.right != null) {
                graph.addEdge(item.x, item.y, item.right.x, item.right.y);
                treeGraphAddEdges(item.right, graph);
            }
        }
    }

    void drawDepthTable(double[] depthTable) {
        if (depthTable != null) {
            System.out.print("\\begin{tabular}{|l");
            for (int i=0; i<depthTable.length/2; i++) {
                System.out.print("|c");
            }
            System.out.println("|}");
            System.out.println("\\hline");
            System.out.print("Old depth:");
            for (int i=0; i<depthTable.length/2; i++) {
                System.out.print("& ");
                System.out.printf("%.3f", depthTable[2*i]);
            }
            System.out.println("\\\\");
            System.out.println("\\hline");
            System.out.print("New depth:");
            for (int i=0; i<depthTable.length/2; i++) {
                System.out.print("& ");
                System.out.printf("%.3f", depthTable[2*i+1]);
            }
            System.out.println("\\\\");
            System.out.println("\\hline");
            System.out.print("\\end{tabular}\n\n");
        }
    }

    void clearDepthList() {
        depthListHead = null;
        depthListEnd = null;
    }

    public static void main (String[] args) {
        Tree tree = new Tree();

        // If no input given, then treat as java Tree 20
        if (args.length == 0) {

            String[] arr = tree.randomPerm(20);
            for (String a : arr) {
                tree.insert(Integer.parseInt(a));
            }

            // y coordinate of root is set to be 100
            tree.findCoordinates(tree.root, 100);

            Graph treeGraph = new Graph();
            tree.treeGraphAddNodes(tree.root, treeGraph);
            tree.treeGraphAddEdges(tree.root, treeGraph);

            treeGraph.outHeader();
            int[] gridCoordinates = treeGraph.getGridCoordinates();
            treeGraph.latex(arr, gridCoordinates, tree);
            treeGraph.outFooter();
            return;
        }

        if (args.length == 1) {
            int n = Integer.parseInt(args[0]);

            if (n < 0) {
                tree.randomTest();
            }

            if (n>0) {
                String[] arr = tree.randomPerm(n);
                for (String a : arr) {
                    tree.insert(Integer.parseInt(a));
                }

                // Assign coordinates to items in tree
                tree.findCoordinates(tree.root, 100);

                // Create graph and add nodes and edges
                Graph treeGraph = new Graph();
                tree.treeGraphAddNodes(tree.root, treeGraph);
                tree.treeGraphAddEdges(tree.root, treeGraph);

                // Output LaTex commands
                treeGraph.outHeader();
                int[] gridCoordinates = treeGraph.getGridCoordinates();
                treeGraph.latex(arr, gridCoordinates, tree);
                treeGraph.outFooter();
            }
            return;
        }

        if (args.length == 2) {

            int n = Integer.parseInt(args[0]);
            int d = Integer.parseInt(args[1]);

            // Inputs need to be positive and d <= n
            if (n < 0 || d < 0) {
                System.out.print("Input positive numbers n & d");
                return;
            }
            if (d > n) {
                d = n;
            }

            String[] perm = tree.randomPerm(n+d);
            Graph treeGraph = new Graph();
            String[] insertedValues = new String[n];

            // First n numbers are inserted
            for (int i=0; i<n; i++) {
                tree.insert(Integer.parseInt(perm[i]));
                insertedValues[i] = perm[i];
            }

            // assign coordinate values to Item class and add to nodeList and edgeList in Graph
            tree.findCoordinates(tree.root, 100);
            tree.treeGraphAddNodes(tree.root, treeGraph);
            tree.treeGraphAddEdges(tree.root, treeGraph);

            treeGraph.outHeader();
            int[] gridCoordinates = treeGraph.getGridCoordinates();
            treeGraph.latex(insertedValues, gridCoordinates, tree);

            // Create new depthList
            tree.clearDepthList();

            // Another random permutation from 1 to n. First d are the indices of the numbers to be deleted
            System.out.println("---- Alternate delete/insert: ---- \n");
            String[] permN = tree.randomPerm(n);

            for (int i=0; i<d; i++) {
                // Delete from tree
                tree.delete(Integer.parseInt(perm[Integer.parseInt(permN[i])]));
                System.out.print("Delete " + perm[Integer.parseInt(permN[i])] + " ");

                // Insert into tree
                tree.insert(Integer.parseInt(perm[i+n]));
                System.out.print("Insert " + perm[i+n] + "\n\n");
            }

            // Create a new graph for a new tree Drawing
            Graph newTreeGraph = new Graph();
            tree.findCoordinates(tree.root, 100);
            tree.treeGraphAddNodes(tree.root, newTreeGraph);
            tree.treeGraphAddEdges(tree.root, newTreeGraph);

            int[] newGridCoordinates = newTreeGraph.getGridCoordinates();
            newTreeGraph.latex(null, newGridCoordinates, tree);
            tree.clearDepthList();

            treeGraph.outFooter();

            return;
        }

        if (args.length == 3) {

            int n = Integer.parseInt(args[0]);
            int d = Integer.parseInt(args[1]);
            int r = Integer.parseInt(args[2]);

            if ( n<0 || d<0 || r<0) {
                System.out.print("Insert positive numbers n d r");
                return;
            }
            if (d > n) {
                d = n;
            }

            // Depth Table to store depths
            double[] depthTable = new double[2*r];

            // treeGraph used to draw tree in first iteration
            Graph treeGraph = new Graph();
            treeGraph.outHeader();

            for (int i=0; i<r; i++) {

                // New tree created for every iteration
                Tree treeRepeat = new Tree();
                String[] perm = treeRepeat.randomPerm(n+d);
                String[] insertedValues = new String[n];

                for(int j=0; j<n; j++) {
                    treeRepeat.insert(Integer.parseInt(perm[j]));
                    insertedValues[j] = perm[j];
                }

                // Create depthList and add average depth to depthTable
                treeRepeat.createDepthList(treeRepeat.root, 0);
                double oldDepth = treeRepeat.averageDepth(treeRepeat.depthListHead, false);
                depthTable[2*i] = oldDepth;
                treeRepeat.clearDepthList();

                // Draw first iteration of tree
                if (i == 0) {
                    treeRepeat.findCoordinates(treeRepeat.root, 100);
                    treeRepeat.treeGraphAddNodes(treeRepeat.root, treeGraph);
                    treeRepeat.treeGraphAddEdges(treeRepeat.root, treeGraph);
                    int[] gridCoordinates = treeGraph.getGridCoordinates();
                    treeGraph.latex(insertedValues, gridCoordinates, treeRepeat);
                    treeRepeat.clearDepthList();
                }

                String[] permN = treeRepeat.randomPerm(n);
                for (int j=0; j<d; j++) {
                    if (i == 0) {
                        System.out.print("Delete " + perm[Integer.parseInt(permN[j])] + " ");
                        System.out.print("Insert " + perm[j+n] + "\n\n");
                    }
                    treeRepeat.delete(Integer.parseInt(perm[Integer.parseInt(permN[j])]));
                    treeRepeat.insert(Integer.parseInt(perm[j+n]));
                }

                // Add average depth of tree to depthTable after insert/delete has happened
                treeRepeat.createDepthList(treeRepeat.root, 0);
                double newDepth = treeRepeat.averageDepth(treeRepeat.depthListHead, false);
                depthTable[2*i+1] = newDepth;
                treeRepeat.clearDepthList();

                // Draw insert/delete tree on first iteration
                if (i == 0) {
                    Graph newTreeGraph = new Graph();
                    treeRepeat.findCoordinates(treeRepeat.root, 100);
                    treeRepeat.treeGraphAddNodes(treeRepeat.root, newTreeGraph);
                    treeRepeat.treeGraphAddEdges(treeRepeat.root, newTreeGraph);

                    int[] newGridCoordinates = newTreeGraph.getGridCoordinates();
                    newTreeGraph.latex(null, newGridCoordinates, treeRepeat);
                    treeRepeat.clearDepthList();
                }

            }

            tree.drawDepthTable(depthTable);
            treeGraph.outFooter();

        }

        if (args.length > 3) {

            for (String arg : args) {
                int x = Integer.parseInt(arg);
                tree.insert(x);
            }

            tree.findCoordinates(tree.root, 100);

            Graph treeGraph = new Graph();
            tree.treeGraphAddNodes(tree.root, treeGraph);
            tree.treeGraphAddEdges(tree.root, treeGraph);

            treeGraph.outHeader();
            int[] gridCoordinates = treeGraph.getGridCoordinates();
            treeGraph.latex(args, gridCoordinates, tree);
            treeGraph.outFooter();

        }
    }
}
