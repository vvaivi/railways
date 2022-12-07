#include "datastructures.hh"

#include <random>

#include <cmath>

#include <set>

#include <queue>

#include <stack>

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

Datastructures::Datastructures()
{
}

Datastructures::~Datastructures()
{
}

/**
 * @brief Datastructures::station_count calculates number of saved stations
 * @return number of stations
 */
unsigned int Datastructures::station_count()
{
    return stations_.size();
}

/**
 * @brief Datastructures::clear_all empties the datastructure containing saved
 * stations
 */
void Datastructures::clear_all()
{
    stations_.clear();
    alphabetically_sorted_stations_.clear();
    sorted_by_distance_stations_.clear();
    station_coordinates_.clear();
    routes_.clear();
}

/**
 * @brief Datastructures::all_stations returns a vector containing ids of saved
 * stations
 * @return station ids
 */
std::vector<StationID> Datastructures::all_stations()
{
    std::vector<StationID> station_ids;
    //avoid reallocating memory
    station_ids.reserve(stations_.size());

    //Saving ids to vector
    std::transform(stations_.begin(), stations_.end(),
                   std::back_inserter(station_ids), [] (const std::unordered_map
                   <StationID,std::shared_ptr<Station>>::value_type &pair)
                   {return pair.first;
    });

    return station_ids;
}

/**
 * @brief Datastructures::add_station saves a new station
 * @param id station id
 * @param name station name
 * @param xy station coordinates
 * @return true if saving successful, false otherwise
 */
bool Datastructures::add_station(StationID id, const Name& name, Coord xy)
{
    //Checking that same station is not saved twice
    if (stations_.find(id) == stations_.end()) {
        //Saving new station
        int distance_to_origo = distance(xy, {0,0});
        std::shared_ptr<Station> new_station(new Station(id, name, xy,
                                                         distance_to_origo));
        stations_.insert({id,new_station});
        station_coordinates_.insert({xy, id});
        return true;
    }

    //Needs to be resorted if something is added
    alphabetically_sorted_stations_.clear();
    sorted_by_distance_stations_.clear();

    //If station with same id is already saved
    return false;
}

/**
 * @brief Datastructures::get_station_name returns the name of given station
 * @param id station id
 * @return station name
 */
Name Datastructures::get_station_name(StationID id)
{
    auto it = stations_.find(id);

    if (it != stations_.end()) {
        return it->second->station_name_;
    }

    //If station with given id does not exist
    return NO_NAME;
}

/**
 * @brief Datastructures::get_station_coordinates returns coordinates of given
 * station
 * @param id station id
 * @return station coordinates
 */
Coord Datastructures::get_station_coordinates(StationID id)
{
    auto it = stations_.find(id);

    if (it != stations_.end()) {
        return it->second->station_coord_;
    }

    //If station with given id does not exist
    return NO_COORD;
}

/**
 * @brief Datastructures::stations_alphabetically returns station ids in
 * alphabetical order
 * @return station ids
 */
std::vector<StationID> Datastructures::stations_alphabetically()
{
    //Checking if order is already saved
    if (alphabetically_sorted_stations_.empty()) {
        std::multimap<std::string,StationID> alphabetical_order;

        //Saving names and ids to multimap
        std::transform(stations_.begin(), stations_.end(), std::inserter(
                       alphabetical_order,alphabetical_order.begin()), [] (const
                       std::pair<StationID,std::shared_ptr<Station>>& pair) {
                       return std::make_pair(pair.second->station_name_,
                                             pair.second->station_id_);
        });

        //Saving ids to vector
        std::transform(alphabetical_order.begin(), alphabetical_order.end(),
                       std::back_inserter(alphabetically_sorted_stations_), []
                       (const std::pair<Name,StationID>& pair) {
                       return pair.second;
        });
    }

    return alphabetically_sorted_stations_;
}

/**
 * @brief Datastructures::stations_distance_increasing returns station ids
 * organized by distance
 * @return station ids
 */
std::vector<StationID> Datastructures::stations_distance_increasing()
{
    //Checking if order is already saved
    if (sorted_by_distance_stations_.empty()) {
        std::multimap<int,StationID> distance_order;

        //Saving names and ids to multimap
        std::transform(stations_.begin(), stations_.end(), std::inserter(
                       distance_order,distance_order.begin()), [] (const
                       std::pair<StationID,std::shared_ptr<Station>>& pair) {
                       return std::make_pair(pair.second->distance_to_origo_,
                                             pair.second->station_id_);
        });

        //Saving ids to vector
        std::transform(distance_order.begin(), distance_order.end(), std::
                       back_inserter(sorted_by_distance_stations_), [] (const
                       std::pair<int,StationID>& pair) {
                       return pair.second;
        });
    }

    return sorted_by_distance_stations_;
}

/**
 * @brief Datastructures::find_station_with_coord returns station id for station
 * in given coord
 * @param xy station coordinates
 * @return station id
 */
StationID Datastructures::find_station_with_coord(Coord xy)
{
    auto it = station_coordinates_.find(xy);

    //If not found
    if (it == station_coordinates_.end()){
        return NO_STATION;
    }

    return it->second;
}

/**
 * @brief Datastructures::change_station_coord changes station coordinates
 * @param id station id
 * @param newcoord station coordinates
 * @return true is successful, false otherwise
 */
bool Datastructures::change_station_coord(StationID id, Coord newcoord)
{
    auto it = stations_.find(id);

    //If not found
    if (it == stations_.end()){
        return false;
    }

    //Updating coordinates
    station_coordinates_.erase(it->second->station_coord_);
    station_coordinates_.insert({newcoord, id});

    //Order and distance changes
    sorted_by_distance_stations_.clear();
    it->second->distance_to_origo_ = distance(newcoord, {0,0});

    it->second->station_coord_ = newcoord;
    return true;
}

/**
 * @brief Datastructures::add_departure saves new departure to a station
 * @param stationid station id
 * @param trainid train id
 * @param time time of departure
 * @return true if successful, false otherwise
 */
bool Datastructures::add_departure(StationID stationid, TrainID trainid, Time
                                   time)
{
    auto station_it = stations_.find(stationid);

    //If station is not found
    if (station_it == stations_.end()) {
        return false;
    }

    //If train with same id and time is already added
    if (std::count(station_it->second->departures_.begin(), station_it->second->
                   departures_.end(), std::make_pair(time, trainid)) == 1) {
        return false;
    }

    station_it->second->departures_.insert(station_it->second->departures_.begin
                                           (), std::make_pair(time, trainid));
    return true;

}

/**
 * @brief Datastructures::remove_departure removes a departure from station
 * @param stationid station id
 * @param trainid train id
 * @param time time of departure
 * @return true if successful, false otherwise
 */
bool Datastructures::remove_departure(StationID stationid, TrainID trainid, Time
                                      time)
{
    auto station_it = stations_.find(stationid);
    auto train_it = find(station_it->second->departures_.begin(),station_it->
                         second->departures_.end(),std::make_pair(time,trainid))
            ;

    //If station or train is not found
    if (station_it == stations_.end() or train_it == station_it->second->
            departures_.end()) {
        return false;
    }

    station_it->second->departures_.erase(train_it);
    return true;
}

/**
 * @brief Datastructures::station_departures_after returns all departures from
 * station after a specific time
 * @param stationid station id
 * @param time time of departure
 * @return true if successful, false otherwise
 */
std::vector<std::pair<Time, TrainID>> Datastructures::station_departures_after
(StationID stationid, Time time)
{
    auto station_it = stations_.find(stationid);

    //If station not found
    if (station_it == stations_.end()) {
        return {{NO_TIME, NO_TRAIN}};
    }

    //Sorting in time order
    std::sort(station_it->second->departures_.begin(),station_it->second->
              departures_.end(), [&] (const std::pair<Time, TrainID> &a, const
              std::pair<Time,TrainID> &b) {
              //If times are equal sorted by id
              if (a.first == b.first) {
                  return a.second < b.second;
              }
              return a.first < b.first;
      });

    //Searching the smallest acceptable time
    auto train_it = std::find_if(station_it->second->departures_.begin(),
                                 station_it->second->departures_.end(), [&]
                                 (const std::pair<Time, TrainID> &train) {
                                 return train.first >= time;
    });

    return {train_it, station_it->second->departures_.end()};
}

/**
 * @brief Datastructures::add_region saves a new region
 * @param id region id
 * @param name region namae
 * @param coords region coordinates
 * @return true if successful, false otherwise
 */
bool Datastructures::add_region(RegionID id, const Name &name,
                                std::vector<Coord> coords)
{
    //Checking that same region is not saved twice
    if (regions_.find(id) == regions_.end()) {
        std::shared_ptr<Region> new_region(new Region(id, name, coords));
        regions_.insert({id,new_region});
        return true;
    }

    //If region with same id is already saved
    return false;
}

/**
 * @brief Datastructures::all_regions returns ids of all saved regions
 * @return region ids
 */
std::vector<RegionID> Datastructures::all_regions()
{
    std::vector<RegionID> region_ids;
    region_ids.reserve(regions_.size());

    //Saving ids to vector
    std::transform(regions_.begin(), regions_.end(),
                   std::back_inserter(region_ids), [] (const std::unordered_map
                   <RegionID,std::shared_ptr<Region>>::value_type &pair)
                   {return pair.first;
    });

    return region_ids;
}

/**
 * @brief Datastructures::get_region_name returns region name
 * @param id region id
 * @return region name
 */
Name Datastructures::get_region_name(RegionID id)
{
    auto it = regions_.find(id);

    if (it != regions_.end()) {
        return it->second->name_;
    }

    //If region with given id does not exist
    return NO_NAME;
}

/**
 * @brief Datastructures::get_region_coords returns coordinates of a region with
 * given id
 * @param id region id
 * @return region coordinates
 */
std::vector<Coord> Datastructures::get_region_coords(RegionID id)
{
    auto it = regions_.find(id);

    if (it != regions_.end()) {
        return it->second->coords_;
    }

    //If region with given id does not exist
    return {NO_COORD};
}

/**
 * @brief Datastructures::add_subregion_to_region saves a new subregion to a
 * region
 * @param id subregion id
 * @param parentid region id
 * @return true if successful, false otherwise
 */
bool Datastructures::add_subregion_to_region(RegionID id, RegionID parentid)
{
    auto it_sub = regions_.find(id);
    auto it_master = regions_.find(parentid);

    //If either region does not exist or region is already a subregion
    if (it_sub == regions_.end() or it_master == regions_.end() or
        it_sub->second->master_region_ != nullptr) {
         return false;
    }

    it_sub->second->master_region_ = &*it_master->second;
    it_master->second->all_subregions_.insert(it_master->second->all_subregions_
                                                .begin(),id);

    //Saving new indirect subregions
    for_each(it_sub->second->all_subregions_.begin(),it_sub->second->
             all_subregions_.end(), [it_master] (RegionID region_id) {
             it_master->second->all_subregions_.insert(it_master->second->
             all_subregions_.begin(),region_id);
    });

    //Subregions need to be added to parents aswell
    Region* parent = it_master->second->master_region_;

    while (parent != nullptr) {
        parent->all_subregions_.insert(parent->all_subregions_.begin(), id);
        parent = parent->master_region_;
    }

    return true;
}

/**
 * @brief Datastructures::add_station_to_region adds station to region
 * @param id station id
 * @param parentid region id
 * @return true if successful, false otherwise
 */
bool Datastructures::add_station_to_region(StationID id, RegionID parentid)
{
    auto it_station = stations_.find(id);
    auto it_region = regions_.find(parentid);

    //If station or region does not exist or station is already in some region
    if (it_station == stations_.end() or it_region == regions_.end() or
        it_station->second->region_ != nullptr) {
         return false;
    }

    it_station->second->region_ = &*it_region->second;
    return true;
}

/**
 * @brief Datastructures::station_in_regions returns all regions that station
 * belongs to
 * @param id station ids
 * @return region ids
 */
std::vector<RegionID> Datastructures::station_in_regions(StationID id)
{
    auto it = stations_.find(id);

    //If station not found
    if (it == stations_.end()) {
        return {NO_REGION};
    }
    //If station is not part of any region
    if (it->second->region_ == nullptr) {
        return {};
    }

    std::vector<RegionID> regions;
    Region* region = it->second->region_;
    regions.insert(regions.end()--,region->region_id_);

    //Looping while master regions exist and saving ids to vector
    while (region->master_region_ != nullptr) {
        regions.insert(regions.end()--,region->master_region_->region_id_);
        region = region->master_region_;
    }

    return regions;
}

/**
 * @brief Datastructures::all_subregions_of_region returns all subregions of
 * a region
 * @param id region id
 * @return subregion ids
 */
std::vector<RegionID> Datastructures::all_subregions_of_region(RegionID id)
{
    auto it = regions_.find(id);

    //If region not found
    if (it == regions_.end()) {
        return {NO_REGION};
    }
    //If region does not have subregion
    if (it->second->all_subregions_.empty()) {
        return {};
    }

    return it->second->all_subregions_;
}

/**
 * @brief Datastructures::stations_closest_to return three closest stations to
 * a given coordinate
 * @param xy coordinate
 * @return station ids
 */
std::vector<StationID> Datastructures::stations_closest_to(Coord xy)
{
    std::multimap<int, StationID> distances;

    //Saving to a multimap to sort
    std::transform(station_coordinates_.begin(), station_coordinates_.end(),
                   std::inserter(distances, distances.begin()), [xy, this]
                   (const std::pair<Coord, StationID>& pair) {
                   return std::make_pair(distance(pair.first, xy),pair.second);
    });

    std::vector<StationID> three_closest;
    three_closest.reserve(3);

    //Saving to a vector
    for (auto& station : distances) {
        if (three_closest.size() == 3) {
            return three_closest;
        }
        three_closest.push_back(station.second);
    }
    return three_closest;
}

/**
 * @brief Datastructures::remove_station removes a station with given id
 * @param id station id
 * @return true if successful, false otherwise
 */
bool Datastructures::remove_station(StationID id)
{
    auto it = stations_.find(id);

    //If station not found
    if(it == stations_.end()) {
        return false;
    }

    //Deleting pointers
    delete it->second->route_back_;
    for (auto& station : it->second->stations_after_) {
        delete station;
    }

    auto it_coord = station_coordinates_.find(it->second->station_coord_);

    //Orders change when something is deleted
    alphabetically_sorted_stations_.clear();
    sorted_by_distance_stations_.clear();

    stations_.erase(it);
    station_coordinates_.erase(it_coord);

    //Deleted from other stations' info
    for (auto& station : it->second->stations_after_) {
        auto it_station = std::find(station->stations_after_.begin(), station->
                                    stations_after_.end(),it->second.get());
        station->stations_after_.erase(it_station);
    }

    //Clearing roads from station
    it->second->stations_after_.clear();

    //Removing from routes
    for (auto& route : routes_){
        auto it = std::find_if(route.second.begin(), route.second.end(), [id]
                               (const std::pair<StationID, Time> & pair) {
                                return pair.first == id;
        });
        if (it != route.second.end()){
            route.second.erase(it);
        }
    }

    //If pointers to station being removed is saved, they are deleted
    for (auto& station : stations_) {
        if (station.second->route_back_->station_id_ == id) {
            delete station.second->route_back_;
        }
        for (auto& station_route : station.second->stations_after_) {
            if (station_route->station_id_ == id) {
                delete station_route;
            }
        }
    }

    return true;
}

/**
 * @brief Datastructures::common_parent_of_regions returns the closest region
 * that both regions belong to
 * @param id1 subregion id
 * @param id2 subregion id
 * @return parent region id
 */
RegionID Datastructures::common_parent_of_regions(RegionID id1, RegionID id2)
{
    auto it_1 = regions_.find(id1);
    auto it_2 = regions_.find(id2);

    //If either region not found
    if (it_1 == regions_.end() or it_2 == regions_.end()) {
        return NO_REGION;
    }

    Region* parent_1 = it_1->second->master_region_;
    Region* parent_2 = it_2->second->master_region_;
    std::vector<RegionID> parents_vector;
    std::set<RegionID> parents_set;

    //Going through parents until same is found
    while (!(parent_1 == nullptr and parent_2 == nullptr)){
        if (parent_1 != nullptr) {
            parents_vector.insert(parents_vector.end()--,parent_1->region_id_);
            parents_set.insert(parent_1->region_id_);
            parent_1 = parent_1->master_region_;
        }
        if (parent_2 != nullptr) {
            parents_vector.insert(parents_vector.end()--,parent_2->region_id_);
            parents_set.insert(parent_2->region_id_);
            parent_2 = parent_2->master_region_;
        }
        if (parents_set.size() != parents_vector.size()) {
            return parents_vector.at(parents_vector.size()-1);
        }
    }

    return NO_REGION;
}

/**
 * @brief Datastructures::distance calculates distance between two points
 * @param coord1 coordinates
 * @param coord2 coordinates
 * @return distance
 */
int Datastructures::distance(Coord coord1, Coord coord2)
{
    return sqrt((coord2.x - coord1.x)*(coord2.x - coord1.x) +
                (coord2.y - coord1.y)*(coord2.y - coord1.y));
}

/**
 * @brief Datastructures::add_train saves a new train
 * @param trainid id
 * @param stationtimes station ids and times
 * @return true if successful, false otherwise
 */
bool Datastructures::add_train(TrainID trainid, std::vector<std::pair<StationID,
                               Time>> stationtimes)
{
    auto it = routes_.find(trainid);

    //If train is already saved
    if (it != routes_.end()) {
        return false;
    }

    for (auto i = stationtimes.begin(); i != stationtimes.end(); i++) {
        //If station not found
        if (stations_.find(i->first) == stations_.end()){
            return false;
        }

        //Saving to departures
        stations_.at(i->first)->departures_.push_back(std::make_pair(i->second,
                                                                     trainid));

        //Saving stations after a particular station on a route and arrival
        //times
        StationID current_station = i->first;
        if (current_station != stationtimes.back().first){
            auto next = std::next(i,1);
            StationID next_station = next->first;

            auto it_current = stations_.find(current_station);
            auto it_next = stations_.find(next_station);

            //Prevent saving twice
            if (std::find(it_current->second->stations_after_.begin(),
                it_current->second->stations_after_.end(),it_next->second.get())
                    == it_current->second->stations_after_.end()){

                it_current->second->stations_after_.push_back(it_next->second.
                                                              get());
            }
        }
    }

    //Saving to routes
    routes_.insert({trainid, stationtimes});

    return true;
}

/**
 * @brief Datastructures::next_stations_from returns ids of stations that come
 * immediately after given station on any route
 * @param id station id
 * @return ids of stations after given station
 */
std::vector<StationID> Datastructures::next_stations_from(StationID id)
{
    auto it = stations_.find(id);

    //If station not found
    if (it == stations_.end()) {
        return {NO_STATION};
    }

    std::vector<StationID> next_stations = {};

    //Saving stations to vector
    for (auto& pair : routes_) {
        auto station_it = std::find_if(pair.second.begin(), pair.second.end(),
                                       [id] (const std::pair<StationID, Time> &
                                       pair) {
                                       return pair.first == id;
        });

        //If next station exists and is not the last station of route
        if (station_it != pair.second.end() and *station_it != pair.second
            .back()) {
            //Moving iterator to next element and saving id
            std::advance(station_it,1);
            next_stations.push_back(station_it->first);
        }
    }

    return next_stations;
}

/**
 * @brief Datastructures::train_stations_from returns all stations that train
 * visits after leaving from given station
 * @param stationid station id
 * @param trainid train id
 * @return ids of stations that train visits
 */
std::vector<StationID> Datastructures::train_stations_from(StationID stationid,
                                                           TrainID trainid)
{
    auto train_it = routes_.find(trainid);

    //If train not found
    if (train_it == routes_.end()) {
        return {NO_STATION};
    }

    auto station_it = std::find_if(train_it->second.begin(), train_it->second.
                                   end(), [stationid] (const std::pair<StationID
                                   , Time> & pair) {
                                   return pair.first == stationid;
    });

    //If station not found
    if (station_it == train_it->second.end()) {
        return {NO_STATION};
    }

    std::vector<StationID> upcoming_stations;

    //Moving iterator to the next element from given station
    std::advance(station_it, 1);

    //Saving station ids to vector
    std::transform(station_it, train_it->second.end(), std::back_inserter(
                   upcoming_stations), [] (const std::pair<StationID, Time> &
                   pair) {
                   return pair.first;
    });

    //If not found
    if (upcoming_stations.empty()) {
        return {NO_STATION};
    }

    return upcoming_stations;
}


/**
 * @brief Datastructures::clear_trains clears all saved trains
 */
void Datastructures::clear_trains()
{
    routes_.clear();

    for (auto& station : stations_) {
        station.second->departures_.clear();
        station.second->stations_after_.clear();
    }
}

/**
 * @brief Datastructures::route_any returns a route between two given stations
 * @param fromid station id of starting station
 * @param toid station id of ending station
 * @return station ids on route and covered distance before station
 */
std::vector<std::pair<StationID, Distance>> Datastructures::route_any(StationID
                                                         fromid, StationID toid)
{
    auto it_from = stations_.find(fromid);
    auto it_to = stations_.find(toid);

    //If either station not found
    if (it_from == stations_.end() or it_to == stations_.end()) {
        return {std::make_pair(NO_STATION, NO_DISTANCE)};
    }

    //Clearing info from last search
    for_each(stations_.begin(), stations_.end(), [] (const std::pair <StationID,
             std::shared_ptr<Station>> & pair) {
             pair.second->colour_ = WHITE;
             pair.second->route_back_ = nullptr;
    });

    //BFS
    it_from->second->colour_ = GREY;

    std::vector<StationID> route;

    std::queue<StationID> help_queue;
    help_queue.push(fromid);

    while(!help_queue.empty()){
        StationID current = help_queue.front();
        help_queue.pop();

        auto it_current = stations_.find(current);

        for(auto& next : it_current->second->stations_after_){
            if(next->colour_ == WHITE){
                next->colour_ = GREY;
                next->route_back_ = it_current->second.get();
                help_queue.push(next->station_id_);
            }
        }

        it_current->second->colour_ = BLACK;
    }

    route.push_back(toid);

    //Searches route from destination to first station
    route = find_path(route, fromid);

    //If route not found
    if(route.size() == 1){
        return {};
    }

    //find_path returns reverse order so needs to be reversed
    std::reverse(route.begin(), route.end());

    //Saving cumulative distances
    std::vector<std::pair<StationID, Distance>> route_with_distances;
    route_with_distances.reserve(route.size());

    int size = route.size();
    route_with_distances.push_back(std::make_pair(route[0], 0));
    for (int i = 1; i < size; i++){
        route_with_distances.push_back(std::make_pair(route[i], distance(
                                       stations_.at(route[i])->station_coord_,
                                       stations_.at(route[i-1])->station_coord_)
                                       + route_with_distances[i-1].second));
    }

    return route_with_distances;
}

/**
 * @brief Datastructures::route_least_stations returns route between two
 * stations that contains least stations
 * @param fromid station id to start from
 * @param toid station id to end to
 * @return station ids on route and covered distance before station
 */
std::vector<std::pair<StationID, Distance>> Datastructures::route_least_stations
(StationID fromid, StationID toid)
{
    //route_any implemented with BFS so returns route that has least stations
    return route_any(fromid, toid);
}

/**
 * @brief Datastructures::route_with_cycle finds a cyclic route
 * @param fromid station id of station that starts and ends the route
 * @return station ids of stations belonging to route
 */
std::vector<StationID> Datastructures::route_with_cycle(StationID fromid)
{
    auto it = stations_.find(fromid);

    //If station not found
    if(it == stations_.end()){
        return {NO_STATION};
    }

    //Clearing info from last search
    for(auto& station : stations_){
        station.second->colour_ = WHITE;
        station.second->route_back_ = nullptr;
    }

    //DFS
    std::stack<StationID> help_stack;
    help_stack.push(fromid);

    std::vector<StationID> route;

    while(!help_stack.empty()){

        StationID current = help_stack.top();
        help_stack.pop();

        auto it_current = stations_.find(current);

        if(it_current->second->colour_ == WHITE){
            it_current->second->colour_ = GREY;
            help_stack.push(it_current->second->station_id_);

            for(const auto& next : it_current->second->stations_after_)
            {
                if (next->colour_ == WHITE){
                    next->route_back_ = it_current->second.get();
                    help_stack.push(next->station_id_);
                }
                //Going through the cycle
                else if(next->colour_ == GREY){
                    //If found in the first element, continued searching
                    if(it_current->second->route_back_ == next){
                        continue;
                    }
                    else {
                        route.push_back(next->station_id_);
                        route.push_back(it_current->second->station_id_);
                        break;
                    }
                }
                if(route.size() == 2){
                    break;
                }
            }
        }
        else {
           it_current->second->colour_ = BLACK;
        }
    }

    //If cycle not found
    if(route.empty()){
        return route;
    }

    route = find_path(route, fromid);

    //Needs to be reversed since find_path returns reversed order
    std::reverse(route.begin(), route.end());

    return route;
}

/**
 * @brief Datastructures::route_shortest_distance finds a route between two
 * station that is the shortest in distance
 * @param fromid station id to start from
 * @param toid station id to end to
 * @return id of stations on route and cumulative route distance
 */
std::vector<std::pair<StationID, Distance>> Datastructures::
route_shortest_distance(StationID fromid, StationID toid)
{
    auto it_from = stations_.find(fromid);
    auto it_to = stations_.find(toid);

    //If either station not found
    if(it_from == stations_.end() or it_to == stations_.end()){
        return {std::make_pair(NO_STATION, NO_DISTANCE)};
    }

    //Clearing info from last search
    for(auto& station : stations_){
        station.second->colour_ = WHITE;
        station.second->route_back_ = nullptr;
        station.second->distance_ = std::numeric_limits<int>::max(); //inf
    }

    //Dijkstran
    it_from->second->colour_ = GREY;
    it_from->second->distance_ = 0;

    std::vector<StationID> route;

    std::set<std::pair<int, StationID>> help_set;
    help_set.insert(std::make_pair(0,fromid));

    while (!help_set.empty()) {
        auto current = help_set.begin();
        help_set.erase(current);

        auto it_current = stations_.find(current->second);

        for(auto& next : it_current->second->stations_after_){
            bool cost_changed = compare_cost_distance(it_current->second.get(),
                                                      next);
            if(next->colour_ == WHITE){
                next->colour_ = GREY;
                next->route_back_ = it_current->second.get();
                help_set.insert(std::make_pair(next->distance_,next
                                               ->station_id_));
            }
            else {
                if (cost_changed) {
                    auto begin = help_set.begin();
                    help_set.erase(begin);
                    help_set.insert(std::make_pair(next->distance_,
                                                   next->station_id_));
                }

            }
        }

        it_current->second->colour_ = BLACK;
    }

    route.push_back(toid);

    //Searches route from destination to first station
    route = find_path(route, fromid);

    //If route not found
    if(route.size() == 1){
        return {};
    }

    //find_path returns reverse order so needs to be reversed
    std::reverse(route.begin(), route.end());

    //Saving cumulative distances
    std::vector<std::pair<StationID, Distance>> route_with_distances;
    route_with_distances.reserve(route.size());

    int size = route.size();
    route_with_distances.push_back(std::make_pair(route[0], 0));
    for (int i = 1; i < size; i++){
        route_with_distances.push_back(std::make_pair(route[i], distance(
                                       stations_.at(route[i])->station_coord_,
                                       stations_.at(route[i-1])->station_coord_)
                                       + route_with_distances[i-1].second));
    }

    return route_with_distances;
}

/**
 * @brief Datastructures::find_path finds path between two routes
 * @param path vector containing station id for first station
 * @param toid station id of last station
 * @return station ids that belong to route
 */
std::vector<StationID> Datastructures::find_path(std::vector<StationID> path,
                                                  StationID toid)
{
    StationID last_station = path.at(path.size()-1);

    auto it = stations_.find(last_station);

    //If not found
    if (it == stations_.end()) {
        return {};
    }

    if(it->second->route_back_ != nullptr){
        path.push_back(it->second->route_back_->station_id_);
        //If destination is found
        if(it->second->route_back_->station_id_ == toid){
            return path;
        }

        return find_path(path, toid);
    }

    return path;
}

/**
 * @brief Datastructures::compare_cost compares and updates costs between routes
 * for route_shortest_distance
 * @param station_a route beginning station
 * @param station_b route ending station
 * @return true if cost is changed, false otherwise
 */
bool Datastructures::compare_cost_distance(Station* station_a,
                                           Station* station_b) {
    int distance_between_stations = distance(station_b->station_coord_,
                                             station_a->station_coord_);

    if (station_b->distance_>station_a->distance_ + distance_between_stations) {
        station_b->distance_ = station_a->distance_ + distance_between_stations;
        station_b->route_back_ = station_a;
        return true;
    }

    return false;
}


