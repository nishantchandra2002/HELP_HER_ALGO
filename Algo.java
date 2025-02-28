import java.util.*;

public class Main {
    static class Point {
        int row;
        int column;
    }

    static double distance(Point p1, Point p2) {
        return Math.sqrt((p1.row - p2.row) * (p1.row - p2.row) + (p1.column - p2.column) * (p1.column - p2.column));
    }

    static double calculatePathDistance(List<Point> path, double[][] graph) {
        double distance = 0.0;
        for (int i = 1; i < path.size(); ++i) {
            distance += graph[path.get(i - 1).row][path.get(i).row];
        }
        return distance;
    }

    static Map.Entry<Double, List<Point>> findShortestPath(double[][] graph, Point source, Point target) {
        int n = graph.length;
        double[] dist = new double[n];
        Arrays.fill(dist, Double.MAX_VALUE);
        int[] prev = new int[n];
        Arrays.fill(prev, -1);
        boolean[] visited = new boolean[n];

        dist[source.row] = 0;

        for (int count = 0; count < n - 1; ++count) {
            int u = -1;
            double minDist = Double.MAX_VALUE;

            for (int i = 0; i < n; ++i) {
                if (!visited[i] && dist[i] < minDist) {
                    u = i;
                    minDist = dist[i];
                }
            }

            if (u == -1 || u == target.row)
                break;

            visited[u] = true;

            for (int v = 0; v < n; ++v) {
                if (!visited[v] && graph[u][v] > 0 && dist[u] + graph[u][v] < dist[v]) {
                    dist[v] = dist[u] + graph[u][v];
                    prev[v] = u;
                }
            }
        }

        List<Point> path = new ArrayList<>();
        int current = target.row;
        while (current != -1) {
            path.add(new Point());
            path.get(path.size() - 1).row = current;
            current = prev[current];
        }
        Collections.reverse(path);

        return new AbstractMap.SimpleEntry<>(calculatePathDistance(path, graph), path);
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        Random rand = new Random();
        System.out.print("Enter the number of rows: ");
        int rows = sc.nextInt();
        System.out.print("Enter the number of columns: ");
        int columns = sc.nextInt();
        double[][] matrix = new double[rows][rows];
        double[][] newMatrix = new double[rows][rows];

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                if (i != j) {
                    matrix[i][j] = 90 + rand.nextInt(100);
                } else {
                    matrix[i][j] = 0.0;
                }
            }
        }

        System.out.println("Distance matrix of points (in meters):");
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                System.out.print(matrix[i][j] + "\t");
                matrix[j][i] = matrix[i][j];
            }
            System.out.println();
        }

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                if (i != j) {
                    newMatrix[i][j] = 90 + rand.nextInt(100);
                } else {
                    newMatrix[i][j] = 0.0;
                }
                newMatrix[j][i] = newMatrix[i][j];
            }
        }

        System.out.println("Traffic Distance Matrix (in meters):");
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                System.out.print(newMatrix[i][j] + "\t");
            }
            System.out.println();
        }

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                matrix[i][j] += newMatrix[i][j];
            }
        }

        System.out.println("Final Distance Matrix (in meters):");
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                System.out.print(matrix[i][j] + "\t");
            }
            System.out.println();
        }

        List<Point> pointA = new ArrayList<>();
        List<Point> pointB = new ArrayList<>();

        for (int i = 0; i < rows / 4; ++i) {
            pointA.add(new Point());
            pointA.get(i).row = rand.nextInt(columns);
            pointA.get(i).column = rand.nextInt(rows);
        }

        for (int i = 0; i < rows / 5; ++i) {
            pointB.add(new Point());
            pointB.get(i).row = rand.nextInt(columns);
            pointB.get(i).column = rand.nextInt(rows);
        }

        System.out.println("Medical Stores:");
        for (Point a : pointA) {
            System.out.println("Point: " + a.row);
        }

        System.out.println("Washrooms:");
        for (Point b : pointB) {
            System.out.println("Point: " + b.row);
        }

        Point currentLocation = new Point();
        System.out.print("Enter your current location (point): ");
        currentLocation.row = sc.nextInt();

        if (currentLocation.row < 0 || currentLocation.row >= rows) {
            System.out.println("Invalid current location.");
            return;
        }

        double shortestDistanceA = Double.MAX_VALUE;
        List<Point> shortestPathA = new ArrayList<>();
        int nearestAIndex = -1;

        for (int i = 0; i < pointA.size(); ++i) {
            double dist = distance(pointA.get(i), currentLocation);
            if (dist < shortestDistanceA) {
                shortestDistanceA = dist;
                nearestAIndex = i;
            }
        }

        if (nearestAIndex != -1) {
            Point source = currentLocation;
            Point target = pointA.get(nearestAIndex);

            Map.Entry<Double, List<Point>> result = findShortestPath(matrix, source, target);
            shortestDistanceA = result.getKey();
            shortestPathA = result.getValue();

            System.out.println("Your current location: Point " + currentLocation.row);
            System.out.println("Nearest Medical Store: Point " + pointA.get(nearestAIndex).row);
            System.out.println("Shortest distance to nearest Medical Store: " + shortestDistanceA + " meters");
            System.out.print("Shortest path to Medical Store: ");
            for (int i = 0; i < shortestPathA.size(); ++i) {
                System.out.print("Point " + shortestPathA.get(i).row);
                if (i < shortestPathA.size() - 1) {
                    System.out.print(" -> ");
                }
            }
            System.out.println();
        } else {
            System.out.println("No nearest Point A found.");
        }

        double shortestDistanceB = Double.MAX_VALUE;
        List<Point> shortestPathB = new ArrayList<>();
        int nearestBIndex = -1;

        for (int i = 0; i < pointB.size(); ++i) {
            double dist = distance(pointB.get(i), pointA.get(nearestAIndex));
            if (dist < shortestDistanceB) {
                shortestDistanceB = dist;
                nearestBIndex = i;
            }
        }

        if (nearestBIndex != -1) {
            Point source = pointA.get(nearestAIndex);
            Point target = pointB.get(nearestBIndex);

            Map.Entry<Double, List<Point>> result = findShortestPath(matrix, source, target);
            shortestDistanceB = result.getKey();
            shortestPathB = result.getValue();

            System.out.println("Nearest Washroom from Medical Store: Point " + pointB.get(nearestBIndex).row);
            System.out.println("Shortest distance from Medical Store to nearest Washroom: " + shortestDistanceB + " meters");
            System.out.print("Shortest path from Medical Store to Washroom: ");
            for (int i = 0; i < shortestPathB.size(); ++i) {
                System.out.print("Point " + shortestPathB.get(i).row);
                if (i < shortestPathB.size() - 1) {
                    System.out.print(" -> ");
                }
            }
            System.out.println();
        } else {
            System.out.println("No nearest Point B from Point A found.");
        }

        System.out.println("Please select a mode of transport");
        System.out.println("Press A for walking");
        System.out.println("Press B for public transport");
        System.out.println("Press C for personal vehicle");
        char x = sc.next().charAt(0);
        double timeToReachPointA = 0.0;
        double timeFromPointAToPointB = 0.0;

        switch (x) {
            case 'A':
                timeToReachPointA = shortestDistanceA / 2.0;
                timeFromPointAToPointB = shortestDistanceB / 2.0;
                System.out.println("Time to reach point A from your current location: " + timeToReachPointA + " seconds");
                System.out.println("Time from point A to point B: " + timeFromPointAToPointB + " seconds");
                break;
            case 'B':
                timeToReachPointA = shortestDistanceA / 15.0;
                timeFromPointAToPointB = shortestDistanceB / 15.0;
                System.out.println("Time to reach point A from your current location: " + timeToReachPointA + " seconds");
                System.out.println("Time from point A to point B: " + timeFromPointAToPointB + " seconds");
                break;
            case 'C':
                timeToReachPointA = shortestDistanceA / 20.0;
                timeFromPointAToPointB = shortestDistanceB / 20.0;
                System.out.println("Time to reach point A from your current location: " + timeToReachPointA + " seconds");
                System.out.println("Time from point A to point B: " + timeFromPointAToPointB + " seconds");
                break;
            default:
                System.out.println("Please choose a mode");
                break;
        }
    }
}
