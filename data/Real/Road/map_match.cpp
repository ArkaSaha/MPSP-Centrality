# include "osrm/match_parameters.hpp"
# include "osrm/nearest_parameters.hpp"
# include "osrm/route_parameters.hpp"
# include "osrm/table_parameters.hpp"
# include "osrm/trip_parameters.hpp"

# include "osrm/coordinate.hpp"
# include "osrm/engine_config.hpp"
# include "osrm/json_container.hpp"

# include "osrm/osrm.hpp"
# include "osrm/status.hpp"

# include <exception>
# include <iostream>
# include <fstream>
# include <filesystem>
# include <string>
# include <utility>
# include <iomanip>

# include <cstdlib>
# include <string>

using namespace std;
using namespace osrm;
using namespace json;

void parse_beijing(char* file, EngineConfig& config)
{
	ifstream trajectory;
	const OSRM osrm{config};
	for (auto& entry : filesystem::recursive_directory_iterator(file))
	{
		string path = entry.path().string();
		if (path.substr(path.size()-4) == ".txt")
		{
			MatchParameters params;
      osrm::engine::api::ResultT result = Object();
			trajectory.open(path);
			string line;
			getline(trajectory,line);
			while (line != "")
			{
				string::size_type p1, p2;
				int id = stoi(line,&p1);
				tm t = {};
				istringstream ss(line.substr(p1+1,19));
				if (ss >> get_time(&t, "%Y-%m-%d %H:%M:%S"))
					params.timestamps.push_back(static_cast<unsigned long>(mktime(&t)));
				double longitude = stod(line.substr(p1+21),&p2);
				double latitude = stod(line.substr(p1+21+p2+1),nullptr);
				params.coordinates.push_back({util::FloatLongitude{longitude},util::FloatLatitude{latitude}});
				getline(trajectory,line);
			}
			trajectory.close();
			params.annotations = true;
			params.annotations_type = RouteParameters::AnnotationsType::All;
			const auto status = osrm.Match(params,result);
			auto& json_result = result.get<Object>();
			if (status == Status::Ok)
			{
				auto& matches = json_result.values["matchings"].get<Array>();
				double max_conf = 0;
				int max_conf_idx = 0;
				for (int j = 0; j < matches.values.size(); j++)
				{
					double conf = matches.values.at(j).get<Object>().values["confidence"].get<Number>().value;
					if (conf > max_conf)
					{
						max_conf = conf;
						max_conf_idx = j;
					}
				}
				for (auto& matching : matches.values.at(max_conf_idx).get<Object>().values["legs"].get<Array>().values)
				{
					auto& match = matching.get<Object>().values["annotation"].get<Object>();
					for (auto& node : match.values["nodes"].get<Array>().values)
						cout << (unsigned long) node.get<Number>().value << "\t";
					cout << endl;
					for (auto& speed : match.values["speed"].get<Array>().values)
						cout << speed.get<Number>().value * 3.6 << "\t";
					cout << endl;
				}
			}
			else
			{
				cerr << path << endl;
				const auto code = json_result.values["code"].get<String>().value;
				const auto message = json_result.values["message"].get<String>().value;
				cerr << "Code: " << code << endl;
				cerr << "Message: " << message << endl;
			}
		}
	}
}

void parse_sf(char* file, EngineConfig& config)
{
	ifstream trajectory;
	const OSRM osrm{config};
	for (auto& entry : filesystem::recursive_directory_iterator(file))
	{
		string path = entry.path().string();
		if (path.substr(path.size()-9) != "_cabs.txt" && path.substr(path.size()-6) != "README")
		{
			cerr << path << endl;
			MatchParameters params;
      osrm::engine::api::ResultT result = Object();
			trajectory.open(path);
			vector< tuple<double,double,unsigned long> > v = vector< tuple<double,double,unsigned long> >();
			string line;
			getline(trajectory,line);
			while (line != "")
			{
				string::size_type p1, p2;
				double latitude = stod(line,&p1);
				double longitude = stod(line.substr(p1+1),&p2);
				params.coordinates.push_back({util::FloatLongitude{longitude},util::FloatLatitude{latitude}});
				unsigned long time = stoul(line.substr(p1+1+p2+3,10),nullptr);
				v.push_back(make_tuple(latitude,longitude,time));
				getline(trajectory,line);
			}
			trajectory.close();
			while (! v.empty())
			{
				double latitude, longitude;
				unsigned long time;
				tie(latitude,longitude,time) = v.back();
				params.coordinates.push_back({util::FloatLongitude{longitude},util::FloatLatitude{latitude}});
				params.timestamps.push_back(time);
				v.pop_back();
			}
			params.annotations = true;
			params.annotations_type = RouteParameters::AnnotationsType::All;
			const auto status = osrm.Match(params,result);
			auto& json_result = result.get<Object>();
			if (status == Status::Ok)
			{
				auto& matches = json_result.values["matchings"].get<Array>();
				double max_conf = 0;
				int max_conf_idx = 0;
				for (int j = 0; j < matches.values.size(); j++)
				{
					double conf = matches.values.at(j).get<Object>().values["confidence"].get<Number>().value;
					if (conf > max_conf)
					{
						max_conf = conf;
						max_conf_idx = j;
					}
				}
				for (auto& matching : matches.values.at(max_conf_idx).get<Object>().values["legs"].get<Array>().values)
				{
					auto& match = matching.get<Object>().values["annotation"].get<Object>();
					for (auto& node : match.values["nodes"].get<Array>().values)
						cout << (unsigned long) node.get<Number>().value << "\t";
					cout << endl;
					for (auto& speed : match.values["speed"].get<Array>().values)
						cout << speed.get<Number>().value * 3.6 << "\t";
					cout << endl;
				}
			}
			else
			{
				cerr << path << endl;
				const auto code = json_result.values["code"].get<String>().value;
				const auto message = json_result.values["message"].get<String>().value;
				cerr << "Code: " << code << endl;
				cerr << "Message: " << message << endl;
			}
		}
	}
}

void parse_porto(char* file, EngineConfig& config)
{
	ifstream trajectory;
	const OSRM osrm{config};
	trajectory.open(strcat(file,"Porto_taxi_data_training.csv"));
	string line;
	getline(trajectory,line);
	getline(trajectory,line);
	while (line != "")
	{
		MatchParameters params;
    osrm::engine::api::ResultT result = Object();
		stringstream lineStream(line);
		string cell;
		int num = 0;
		vector<string> res = vector<string>();
		while (num < 8 && getline(lineStream,cell,','))
		{
			num++;
			res.push_back(cell);
		}
		unsigned long time = stoul(res[5],nullptr);
		params.timestamps.push_back(time);
		while (getline(lineStream,cell,']') && cell != "" && cell != "[")
		{
			string::size_type p;
			double longitude = stod(cell.substr(3),&p);
			double latitude = stod(cell.substr(p+5),nullptr);
			params.coordinates.push_back({util::FloatLongitude{longitude},util::FloatLatitude{latitude}});
			params.timestamps.push_back(time);
			time += 15;
		}
		params.annotations = true;
		params.annotations_type = RouteParameters::AnnotationsType::All;
		const auto status = osrm.Match(params,result);
		auto& json_result = result.get<Object>();
		if (status == Status::Ok)
		{
			auto& matches = json_result.values["matchings"].get<Array>();
			double max_conf = 0;
			int max_conf_idx = 0;
			for (int j = 0; j < matches.values.size(); j++)
			{
				double conf = matches.values.at(j).get<Object>().values["confidence"].get<Number>().value;
				if (conf > max_conf)
				{
					max_conf = conf;
					max_conf_idx = j;
				}
			}
			for (auto& matching : matches.values.at(max_conf_idx).get<Object>().values["legs"].get<Array>().values)
			{
				auto& match = matching.get<Object>().values["annotation"].get<Object>();
				for (auto& node : match.values["nodes"].get<Array>().values)
					cout << (unsigned long) node.get<Number>().value << "\t";
				cout << endl;
				for (auto& speed : match.values["speed"].get<Array>().values)
					cout << speed.get<Number>().value * 3.6 << "\t";
				cout << endl;
			}
		}
		else
		{
			const auto code = json_result.values["code"].get<String>().value;
			const auto message = json_result.values["message"].get<String>().value;
			cerr << "Code: " << code << endl;
			cerr << "Message: " << message << endl;
		}
		getline(trajectory,line);
	}
	trajectory.close();
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cerr << "Usage : " << argv[0] << " data.osrm path-to-trajectories" << endl;
		return EXIT_FAILURE;
	}
	EngineConfig config;
	config.storage_config = {argv[1]};
	config.use_shared_memory = false;
	config.algorithm = EngineConfig::Algorithm::MLD;
	parse_porto(argv[2],config);
	// parse_beijing(argv[2],config);
	// parse_sf(argv[2],config);
	return EXIT_SUCCESS;
}
