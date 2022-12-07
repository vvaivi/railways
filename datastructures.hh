#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH

#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <functional>
#include <exception>
#include <memory>
#include <iostream>
#include <map>
#include <unordered_set>

// Types for IDs
using StationID = std::string;
using TrainID = std::string;
using RegionID = unsigned long long int;
using Name = std::string;
using Time = unsigned short int;

// Return values for cases where required thing was not found
StationID const NO_STATION = "---";
TrainID const NO_TRAIN = "---";
RegionID const NO_REGION = -1;
Name const NO_NAME = "!NO_NAME!";
Time const NO_TIME = 9999;

// Return value for cases where integer values were not found
int const NO_VALUE = std::numeric_limits<int>::min();

// Type for a coordinate (x, y)
struct Coord
{
    int x = NO_VALUE;
    int y = NO_VALUE;
};

// Example: Defining == and hash function for Coord so that it can be used
// as key for std::unordered_map/set, if needed
inline bool operator==(Coord c1, Coord c2) { return c1.x == c2.x && c1.y == c2.y; }
inline bool operator!=(Coord c1, Coord c2) { return !(c1==c2); } // Not strictly necessary

struct CoordHash
{
    std::size_t operator()(Coord xy) const
    {
        auto hasher = std::hash<int>();
        auto xhash = hasher(xy.x);
        auto yhash = hasher(xy.y);
        // Combine hash values (magic!)
        return xhash ^ (yhash + 0x9e3779b9 + (xhash << 6) + (xhash >> 2));
    }
};

// Example: Defining < for Coord so that it can be used
// as key for std::map/set
inline bool operator<(Coord c1, Coord c2)
{
    if (c1.y < c2.y) { return true; }
    else if (c2.y < c1.y) { return false; }
    else { return c1.x < c2.x; }
}

// Return value for cases where coordinates were not found
Coord const NO_COORD = {NO_VALUE, NO_VALUE};

// Type for a distance (in metres)
using Distance = int;

// Return value for cases where Distance is unknown
Distance const NO_DISTANCE = NO_VALUE;

// This exception class is there just so that the user interface can notify
// about operations which are not (yet) implemented
class NotImplemented : public std::exception
{
public:
    NotImplemented() : msg_{} {}
    explicit NotImplemented(std::string const& msg) : msg_{msg + " not implemented"} {}

    virtual const char* what() const noexcept override
    {
        return msg_.c_str();
    }
private:
    std::string msg_;
};

// This is the class you are supposed to implement

class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    // Estimate of performance: O(1)
    // Short rationale for estimate: STL documentation for size (constant)
    unsigned int station_count();

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for clear (linear)
    void clear_all();

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for transform (linear)
    std::vector<StationID> all_stations();

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear), clear
    // (linear) and insert (linear)
    bool add_station(StationID id, Name const& name, Coord xy);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    Name get_station_name(StationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    Coord get_station_coordinates(StationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for transform (linear)
    std::vector<StationID> stations_alphabetically();

    // Estimate of performance: O(nlog(n))
    // Short rationale for estimate: STL documentation for sort (linearithmic)
    std::vector<StationID> stations_distance_increasing();

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    StationID find_station_with_coord(Coord xy);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find and clear
    // (linear)
    bool change_station_coord(StationID id, Coord newcoord);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find and count
    // (linear)
    bool add_departure(StationID stationid, TrainID trainid, Time time);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find and erase
    // (linear)
    bool remove_departure(StationID stationid, TrainID trainid, Time time);

    // Estimate of performance: O(nlog(n))
    // Short rationale for estimate: STL documentation for sort (linearithmic)
    std::vector<std::pair<Time, TrainID>> station_departures_after(StationID stationid, Time time);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    bool add_region(RegionID id, Name const& name, std::vector<Coord> coords);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for transform (linear)
    std::vector<RegionID> all_regions();

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    Name get_region_name(RegionID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    std::vector<Coord> get_region_coords(RegionID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear) and
    // for and while loops
    bool add_subregion_to_region(RegionID id, RegionID parentid);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    bool add_station_to_region(StationID id, RegionID parentid);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear) and
    // while loop
    std::vector<RegionID> station_in_regions(StationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear)
    std::vector<RegionID> all_subregions_of_region(RegionID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for transform (linear)
    // and for loop
    std::vector<StationID> stations_closest_to(Coord xy);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: STL documentation for find and erase
    // (linear) used inside a for loop
    bool remove_station(StationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL documentation for find (linear) and
    // while loop
    RegionID common_parent_of_regions(RegionID id1, RegionID id2);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: STL-documentation for find and find_if
    // (linear), inside a for loop makes n²
    bool add_train(TrainID trainid, std::vector<std::pair<StationID, Time>> stationtimes);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: STL-documentation for find_if (linear)
    // used inside a for loop
    std::vector<StationID> next_stations_from(StationID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: STL-documentation for transform, find and
    // find_if (linear)
    std::vector<StationID> train_stations_from(StationID stationid, TrainID trainid);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: STL-documentation for clear (linear) used
    // inside for loop
    void clear_trains();

    // Estimate of performance: O(n²)
    // Short rationale for estimate: For loop inside while loop and find_path
    // which is also O(n²)
    std::vector<std::pair<StationID, Distance>> route_any(StationID fromid, StationID toid);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: For loop inside while loop and find_path
    // which is also O(n²)
    std::vector<std::pair<StationID, Distance>> route_least_stations(StationID fromid, StationID toid);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: For loop inside while loop and  find_path
    // which is also O(n²)
    std::vector<StationID> route_with_cycle(StationID fromid);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: For loop inside while loop and  find_path
    // which is also O(n²)
    std::vector<std::pair<StationID, Distance>> route_shortest_distance(StationID fromid, StationID toid);

private:
    // For BFS and DFS and Dijkstran
    enum Colour {WHITE, GREY, BLACK};

    struct Region {
        RegionID region_id_;
        Name name_;
        std::vector<Coord> coords_;

        Region* master_region_ = nullptr;

        //Direct and indirect subregions
        std::vector<RegionID> all_subregions_;

        Region(RegionID id, Name name, std::vector<Coord> coords){
            region_id_ = id;
            name_ = name;
            coords_ = coords;
        }
    };

    struct Station {
        StationID station_id_;
        Name station_name_;
        Coord station_coord_;
        int distance_to_origo_;

        Region* region_ = nullptr;
        std::vector<std::pair<Time, TrainID>> departures_;

        //Stations that come after this station on some route
        std::vector<Station*> stations_after_;

        //For BFS and DFS and Dijkstran
        Colour colour_;
        Station* route_back_;

        //For Dijkstran
        int distance_;

        Station (StationID id, Name name, Coord coord, int distance) {
            station_id_ = id;
            station_name_ = name;
            station_coord_ = coord;
            distance_to_origo_ = distance;
        }
    };

    std::unordered_map<StationID,std::shared_ptr<Station>> stations_;
    std::unordered_map<RegionID,std::shared_ptr<Region>> regions_;
    std::unordered_map<TrainID, std::vector<std::pair<StationID, Time>>>
    routes_;

    //To avoid sorting every time when method called
    std::vector<StationID> alphabetically_sorted_stations_;
    std::vector<StationID> sorted_by_distance_stations_;

    //Saving coordinates and ids separately to make methods more efficient
    std::unordered_map<Coord, StationID, CoordHash> station_coordinates_;

    // Calculates the distance between two points
    // Estimate of performance: O(1)
    // Short rationale for estimate: only calculation
    int distance(Coord coord1, Coord coord2);

    // Estimate of performance: O(n²)
    // Short rationale for estimate: find (linear) and recursive function
    std::vector<StationID> find_path(std::vector<StationID> path, StationID toid
                                     );

    // Estimate of performance: O(1)
    // Short rationale for estimate: only comparison and calculation
    bool compare_cost_distance(Station* station_a, Station* station_b);
};

#endif // DATASTRUCTURES_HH
