#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>
using namespace std;

struct Point {
    int row;
    int column;
};

double distance(const Point& p1, const Point& p2) {
    return sqrt((p1.row - p2.row) * (p1.row - p2.row) + (p1.column - p2.column) * (p1.column - p2.column));
}

double calculatePathDistance(const vector<Point>& path, const vector<vector<double>>& graph) {
    double distance = 0.0;
    for (int i = 1; i < path.size(); ++i) {
        distance += graph[path[i - 1].row][path[i].row];
    }
    return distance;
}

pair<double, vector<Point>> findShortestPath(
    const vector<vector<double>>& graph, const Point& source, const Point& target) {

    int n = graph.size();
    vector<double> dist(n, numeric_limits<double>::max());
    vector<int> prev(n, -1);
    vector<bool> visited(n, false);

    dist[source.row] = 0;

    for (int count = 0; count < n - 1; ++count) {
        int u = -1;
        double minDist = numeric_limits<double>::max();

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
    vector<Point> path;
    int current = target.row;
    while (current != -1) {
        path.push_back(Point{current, 0}); // In single point form, the column is not needed.
        current = prev[current];
    }
    reverse(path.begin(), path.end());

    return make_pair(calculatePathDistance(path, graph), path);
}

int main() {
    srand(static_cast<unsigned>(time(0)));
    int rows, columns;
    cout << "Enter the number of rows: ";
    cin >> rows;
    cout << "Enter the number of columns: ";
    cin >> columns;
    vector<vector<double>> matrix(rows, vector<double>(rows));
    vector<vector<double>> newMatrix(rows, vector<double>(rows));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            if (i != j) {
                matrix[i][j] = 90 + static_cast<double>(rand() % 100); // Ensure distances are greater than 90 meters
            } else {
                matrix[i][j] = 0.0;
            }
        }
    }
    
    cout << "Distance matrix of points (in meters):" << endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            cout << matrix[i][j] << "\t";
            matrix[j][i] = matrix[i][j];
        }
        
        cout << endl;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            if (i != j) {
                newMatrix[i][j] = 90 + static_cast<double>(rand() % 100); // Ensure distances are greater than 90 meters
            } else {
                newMatrix[i][j] = 0.0; 
            }
            newMatrix[j][i] = newMatrix[i][j]; 
        }
    }
    
    cout << "Traffic Distance Matrix (in meters):" << endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            cout << newMatrix[i][j] << "\t";
        }
        cout << endl;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            matrix[i][j] += newMatrix[i][j];
        }
    }
    
    cout << "Final Distance Matrix (in meters):" << endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
    vector<Point> pointA;
    vector<Point> pointB;

    for (int i = 0; i < rows / 4; ++i) {
        pointA.push_back(Point{rand() % columns, rand() % rows});
    }

    for (int i = 0; i < rows / 5; ++i) {
        pointB.push_back(Point{rand() % columns, rand() % rows});
    }

    cout << "Medical Stores:" << endl;
    for (const auto& a : pointA) {
        cout << "Point: " << a.row << endl;
    }

    cout << "Washrooms:" << endl;
    for (const auto& b : pointB) {
        cout << "Point: " << b.row << endl;
    }
    Point currentLocation;
    cout << "Enter your current location (point): ";
    cin >> currentLocation.row;

    if (currentLocation.row < 0 || currentLocation.row >= rows) {
        cout << "Invalid current location." << endl;
        return 1;
    }
    double shortestDistanceA = numeric_limits<double>::max();
    vector<Point> shortestPathA;
    int nearestAIndex = -1;

    for (int i = 0; i < pointA.size(); ++i) {
        double dist = distance(pointA[i], currentLocation);
        if (dist < shortestDistanceA) {
            shortestDistanceA = dist;
            nearestAIndex = i;
        }
    }

    if (nearestAIndex != -1) {
        Point source = currentLocation;
        Point target = pointA[nearestAIndex];

        pair<double, vector<Point>> result = findShortestPath(matrix, source, target);
        shortestDistanceA = result.first;
        shortestPathA = result.second;

        cout << "Your current location: Point " << currentLocation.row << endl;
        cout << "Nearest Medical Store: Point " << pointA[nearestAIndex].row << endl;
        cout << "Shortest distance to nearest Medical Store: " << shortestDistanceA << " meters" << endl;
        cout << "Shortest path to Medical Store: ";
        for (int i = 0; i < shortestPathA.size(); ++i) {
            cout << "Point " << shortestPathA[i].row;
            if (i < shortestPathA.size() - 1) {
                cout << " -> ";
            }
        }
        cout << endl;
    } else {
        cout << "No nearest Point A found." << endl;
    }
    double shortestDistanceB = numeric_limits<double>::max();
    vector<Point> shortestPathB;
    int nearestBIndex = -1;

    for (int i = 0; i < pointB.size(); ++i) {
        double dist = distance(pointB[i], pointA[nearestAIndex]);
        if (dist < shortestDistanceB) {
            shortestDistanceB = dist;
            nearestBIndex = i;
        }
    }

    if (nearestBIndex != -1) {
        Point source = pointA[nearestAIndex];
        Point target = pointB[nearestBIndex];

        pair<double, vector<Point>> result = findShortestPath(matrix, source, target);
        shortestDistanceB = result.first;
        shortestPathB = result.second;

        cout << "Nearest Washroom from Medical Store: Point " << pointB[nearestBIndex].row << endl;
        cout << "Shortest distance from Medical Store to nearest Washroom: " << shortestDistanceB << " meters" << endl;
        cout << "Shortest path from Medical Store to Washroom: ";
        for (int i = 0; i < shortestPathB.size(); ++i) {
            cout << "Point " << shortestPathB[i].row;
            if (i < shortestPathB.size() - 1) {
                cout << " -> ";
            }
        }
        cout << endl;
    } else {
        cout << "No nearest Point B from Point A found." << endl;
    }
    char x;
    cout << "Please select a mode of transport" << endl << "Press A for walking" << endl << "Press B for public transport" << endl << "Press C for personal vehicle" << endl;
    cin >> x;
    double timeToReachPointA = 0.0;
    double timeFromPointAToPointB = 0.0;
    
    switch (x) {
    case 'A':
        timeToReachPointA = shortestDistanceA / 2.0;
        timeFromPointAToPointB = shortestDistanceB / 2.0;
        cout << "Time to reach point A from your current location: " << timeToReachPointA << " seconds" << endl;
        cout << "Time from point A to point B: " << timeFromPointAToPointB << " seconds" << endl;
        break;
    case 'B':
    timeToReachPointA = shortestDistanceA /15.0;
     timeFromPointAToPointB = shortestDistanceB / 15.0;
    cout << "Time to reach point A from your current location: " << timeToReachPointA << " seconds" << endl;
    cout << "Time from point A to point B: " << timeFromPointAToPointB << " seconds" << endl;
    break;
    case 'C':
    timeToReachPointA = shortestDistanceA / 20.0;
    timeFromPointAToPointB = shortestDistanceB / 20.0;
    cout << "Time to reach point A from your current location: " << timeToReachPointA << " seconds" << endl;
    cout << "Time from point A to point B: " << timeFromPointAToPointB << " seconds" << endl;
    break;
    default: 
        cout << "Please choose a mode"; 
        break;
    }
    return 0;
}
