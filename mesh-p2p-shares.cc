/*
* Mesh/Random P2P Network with Share Distribution
* Created by: Assistant for malek-attia
* Date: 2025-04-19
*/

#include "ns3/applications-module.h"
#include "ns3/core-module.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/netanim-module.h"
#include "ns3/network-module.h"
#include "ns3/point-to-point-module.h"

#include <algorithm>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("MeshP2PSharesExample");

// SharesApp Class Definition
class SharesApp : public Application
{
public:
    // Proper TypeId registration - note we're using "ns3::SharesApp" here
    static TypeId GetTypeId()
    {
        static TypeId tid = TypeId("ns3::SharesApp")
                                .SetParent<Application>()
                                .SetGroupName("P2PShares")
                                .AddConstructor<SharesApp>()
                                // Add trace source for packet transmission
                                .AddTraceSource("Tx", 
                                            "A packet has been sent",
                                            MakeTraceSourceAccessor(&SharesApp::m_txTrace),
                                            "ns3::Packet::TracedCallback");
        return tid;
    }

    SharesApp()
        : m_sendEvent(),
        m_socket(nullptr),
        m_packetsSent(0),
        m_packetSize(1024),
        m_sendInterval(Seconds(1.0)),
        m_runningTime(Seconds(10.0)),
        m_nodeId(0),
        m_shareCount(0)
    {
    }

    virtual ~SharesApp()
    {
        m_socket = nullptr;
    }

    void Setup(uint32_t nodeId,
                std::vector<Ipv4Address> peerAddresses,
                std::vector<uint16_t> peerPorts)
    {
        m_nodeId = nodeId;
        m_peerAddresses = peerAddresses;
        m_peerPorts = peerPorts;
    }

    void SetSendInterval(Time interval)
    {
        m_sendInterval = interval;
    }

    void SetPacketSize(uint32_t size)
    {
        m_packetSize = size;
    }

    void SetRunningTime(Time time)
    {
        m_runningTime = time;
    }

protected:
    virtual void DoDispose()
    {
        m_socket = nullptr;
        Application::DoDispose();
    }

private:
    virtual void StartApplication()
    {
        // Create the socket with proper binding
        m_socket = Socket::CreateSocket(GetNode(), UdpSocketFactory::GetTypeId());
        InetSocketAddress local = InetSocketAddress(Ipv4Address::GetAny(), 0); // Use any available port
        m_socket->Bind(local);

        NS_LOG_INFO("Node " << m_nodeId << " sender started with " 
                << m_peerAddresses.size() << " peers");

        // Schedule the first send event
        m_sendEvent = Simulator::Schedule(Seconds(1.0), &SharesApp::SendShare, this);

        // Schedule the stop event
        Simulator::Schedule(m_runningTime, &SharesApp::StopApplication, this);
    }

    virtual void StopApplication()
    {
        Simulator::Cancel(m_sendEvent);
        if (m_socket)
        {
            m_socket->Close();
        }
        
        NS_LOG_INFO("Node " << m_nodeId << " sender stopped after sending " 
                << m_packetsSent << " packets");
    }

    void SendShare()
    {
        NS_ASSERT(m_sendEvent.IsExpired());

        if (m_peerAddresses.empty())
        {
            NS_LOG_INFO("Node " << m_nodeId << " has no peers to send shares to.");
            return;
        }

        // Add debug output
        NS_LOG_INFO("Node " << m_nodeId << " attempting to send share at time " 
                    << Simulator::Now().GetSeconds() << "s");
        
        // Select a random peer to send the share to
        Ptr<UniformRandomVariable> random = CreateObject<UniformRandomVariable>();
        uint32_t peerIndex = random->GetInteger(0, m_peerAddresses.size() - 1);
        Ipv4Address peerAddress = m_peerAddresses[peerIndex];
        uint16_t peerPort = m_peerPorts[peerIndex];

        // Create share message
        std::stringstream shareMsg;
        shareMsg << "SHARE_" << m_nodeId << "_" << m_shareCount << "_"
                << Simulator::Now().GetSeconds();
        std::string shareMsgStr = shareMsg.str();

        Ptr<Packet> packet = Create<Packet>((uint8_t*)shareMsgStr.c_str(), shareMsgStr.length());

        // Send the share
        int result = m_socket->SendTo(packet, 0, InetSocketAddress(peerAddress, peerPort));
        
        // Add debug output for send result
        if (result >= 0) {
            NS_LOG_INFO("Node " << m_nodeId << " successfully sent share #" << m_shareCount 
                    << " to peer at " << peerAddress << ":" << peerPort);
            m_packetsSent++;
            m_shareCount++;
            
            // Trigger the trace source for packet transmission
            m_txTrace(packet);
        } else {
            NS_LOG_INFO("Node " << m_nodeId << " failed to send share to peer at " 
                    << peerAddress << ":" << peerPort << " with error code " << result);
        }

        // Schedule the next share to be sent
        m_sendEvent = Simulator::Schedule(m_sendInterval, &SharesApp::SendShare, this);
    }

    EventId m_sendEvent;
    Ptr<Socket> m_socket;
    uint32_t m_packetsSent;
    uint32_t m_packetSize;
    Time m_sendInterval;
    Time m_runningTime;
    uint32_t m_nodeId;
    std::vector<Ipv4Address> m_peerAddresses;
    std::vector<uint16_t> m_peerPorts;
    uint32_t m_shareCount;
    
    // Trace source for packet transmission
    TracedCallback<Ptr<const Packet>> m_txTrace;
};

// SharesReceiver Class Definition
class SharesReceiver : public Application
{
public:
    // Proper TypeId registration - note we're using "ns3::SharesReceiver" here
    static TypeId GetTypeId()
    {
        static TypeId tid = TypeId("ns3::SharesReceiver")
                                .SetParent<Application>()
                                .SetGroupName("P2PShares")
                                .AddConstructor<SharesReceiver>();
        return tid;
    }

    SharesReceiver()
        : m_socket(nullptr),
        m_port(0),
        m_packetsReceived(0),
        m_nodeId(0)
    {
    }

    virtual ~SharesReceiver()
    {
        m_socket = nullptr;
    }

    void Setup(uint16_t port, uint32_t nodeId)
    {
        m_port = port;
        m_nodeId = nodeId;
    }

protected:
    virtual void DoDispose()
    {
        m_socket = nullptr;
        Application::DoDispose();
    }

private:
    virtual void StartApplication()
    {
        // Create the socket
        m_socket = Socket::CreateSocket(GetNode(), UdpSocketFactory::GetTypeId());
        InetSocketAddress local = InetSocketAddress(Ipv4Address::GetAny(), m_port);
        m_socket->Bind(local);
        m_socket->SetRecvCallback(MakeCallback(&SharesReceiver::HandleRead, this));
        
        NS_LOG_INFO("Node " << m_nodeId << " receiver started at port " << m_port);
    }

    virtual void StopApplication()
    {
        if (m_socket)
        {
            m_socket->Close();
            m_socket->SetRecvCallback(MakeNullCallback<void, Ptr<Socket>>());
        }
        
        NS_LOG_INFO("Node " << m_nodeId << " receiver stopped after receiving " 
                << m_packetsReceived << " packets");
    }

    void HandleRead(Ptr<Socket> socket)
    {
        Ptr<Packet> packet;
        Address from;
        while ((packet = socket->RecvFrom(from)))
        {
            m_packetsReceived++;

            // Extract the share message from the packet
            uint8_t buf[packet->GetSize() + 1];
            packet->CopyData(buf, packet->GetSize());
            buf[packet->GetSize()] = '\0';
            std::string message = std::string((char*)buf);

            // Get sender information
            InetSocketAddress senderAddress = InetSocketAddress::ConvertFrom(from);
            Ipv4Address senderIp = senderAddress.GetIpv4();
            uint16_t senderPort = senderAddress.GetPort();

            NS_LOG_INFO("Node " << m_nodeId << " received share: \"" << message << "\" from "
                        << senderIp << ":" << senderPort << " at time "
                        << Simulator::Now().GetSeconds() << "s");
        }
    }

    Ptr<Socket> m_socket;
    uint16_t m_port;
    uint32_t m_packetsReceived;
    uint32_t m_nodeId;
};

// Simplified connection tracking to store IP addresses directly during creation
struct NodeConnection {
    uint32_t node1;
    uint32_t node2;
    Ipv4Address addr1;  // Address of node1 in this connection
    Ipv4Address addr2;  // Address of node2 in this connection
};

// Packet Tx Callback function
static void
PacketTxCallback(Ptr<const Packet> packet)
{
    NS_LOG_INFO("PACKET TX: " << packet->GetSize() << " bytes at time " 
            << Simulator::Now().GetSeconds() << "s");
}

int
main(int argc, char* argv[])
{
    // Set the random seed for reproducibility
    RngSeedManager::SetSeed(1);
    
    // DO NOT set custom scheduler - use default
    // Simulator::SetScheduler(ObjectFactory("ns3::MapScheduler"));

    // Configuration parameters
    uint32_t nNodes = 12;       // Number of nodes in the network
    uint32_t nPeers = 3;        // Number of peers per node
    double sendInterval = 2.0;  // Interval between sending shares (seconds)
    double simTime = 30.0;      // Simulation time (seconds)
    double nodeSpacing = 100.0; // Distance between nodes for visualization
    uint16_t basePort = 9000;   // Base port number for applications

    // Command line arguments
    CommandLine cmd(__FILE__);
    cmd.AddValue("nNodes", "Number of nodes", nNodes);
    cmd.AddValue("nPeers", "Number of peers per node", nPeers);
    cmd.AddValue("sendInterval", "Interval between sending shares (seconds)", sendInterval);
    cmd.AddValue("simTime", "Simulation time (seconds)", simTime);
    cmd.Parse(argc, argv);

    // Validate parameters
    if (nPeers >= nNodes)
    {
        nPeers = nNodes - 1;
        NS_LOG_WARN("Number of peers adjusted to " << nPeers << " (maximum possible)");
    }

    // Enable logging
    Time::SetResolution(Time::NS);
    LogComponentEnable("MeshP2PSharesExample", LOG_LEVEL_INFO);
    LogComponentEnable("UdpSocketImpl", LOG_LEVEL_INFO);

    // Create nodes
    NS_LOG_INFO("Creating " << nNodes << " nodes");
    NodeContainer nodes;
    nodes.Create(nNodes);

    // Install internet stack on all nodes
    InternetStackHelper internet;
    internet.Install(nodes);

    // Define point-to-point link properties with hard-coded latency
    PointToPointHelper p2p;
    p2p.SetDeviceAttribute("DataRate", StringValue("5Mbps"));
    p2p.SetChannelAttribute("Delay", StringValue("10ms")); // Hard-coded latency

    // Connection tracking to easily map nodes to IP addresses
    std::vector<NodeConnection> connections;
    std::vector<std::vector<uint32_t>> peerConnections(nNodes);

    // Random number generator for peer selection
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a connection topology
    NS_LOG_INFO("Creating random P2P topology");

    // Count connections per node
    std::vector<uint32_t> connectionCount(nNodes, 0);

    // Create the initial topology
    for (uint32_t i = 0; i < nNodes; i++)
    {
        // Create a list of potential peers (all other nodes)
        std::vector<uint32_t> potentialPeers;
        for (uint32_t j = 0; j < nNodes; j++)
        {
            if (i != j)
            {
                potentialPeers.push_back(j);
            }
        }

        // Shuffle the list of potential peers
        std::shuffle(potentialPeers.begin(), potentialPeers.end(), gen);

        // Connect to a number of peers
        uint32_t peersToConnect = std::min(nPeers, static_cast<uint32_t>(potentialPeers.size()));

        for (uint32_t p = 0; p < peersToConnect; p++)
        {
            uint32_t peer = potentialPeers[p];

            // Check if this connection already exists (to avoid duplicates)
            bool connectionExists = false;
            for (uint32_t c = 0; c < peerConnections[peer].size(); c++)
            {
                if (peerConnections[peer][c] == i)
                {
                    connectionExists = true;
                    break;
                }
            }

            if (!connectionExists)
            {
                peerConnections[i].push_back(peer);
                connectionCount[i]++;
                connectionCount[peer]++;

                // Create a new NodeContainer for this p2p link
                NodeContainer nc;
                nc.Add(nodes.Get(i));
                nc.Add(nodes.Get(peer));

                // Install the p2p link between these nodes
                NetDeviceContainer ndc = p2p.Install(nc);

                // Assign IP addresses for this link
                Ipv4AddressHelper ipv4;
                std::ostringstream subnet;
                subnet << "10." << (i + 1) << "." << (peer + 1) << ".0";
                ipv4.SetBase(subnet.str().c_str(), "255.255.255.0");
                Ipv4InterfaceContainer ifc = ipv4.Assign(ndc);

                // Track connections with their addresses
                NodeConnection conn;
                conn.node1 = i;
                conn.node2 = peer;
                conn.addr1 = ifc.GetAddress(0);
                conn.addr2 = ifc.GetAddress(1);
                connections.push_back(conn);

                NS_LOG_INFO("Created connection between Node " << i << " and Node " << peer 
                        << " with addresses " << conn.addr1 << " and " << conn.addr2);
            }
        }
    }

    // Ensure all nodes have at least one connection
    for (uint32_t i = 0; i < nNodes; i++)
    {
        if (connectionCount[i] == 0)
        {
            // Find a random node to connect to
            std::uniform_int_distribution<> dis(0, nNodes - 1);
            uint32_t peer;
            do
            {
                peer = dis(gen);
            } while (peer == i);

            // Create a new NodeContainer for this p2p link
            NodeContainer nc;
            nc.Add(nodes.Get(i));
            nc.Add(nodes.Get(peer));

            // Install the p2p link between these nodes
            NetDeviceContainer ndc = p2p.Install(nc);

            // Assign IP addresses for this link
            Ipv4AddressHelper ipv4;
            std::ostringstream subnet;
            subnet << "10." << (i + 1) << "." << (peer + 1) << ".0";
            ipv4.SetBase(subnet.str().c_str(), "255.255.255.0");
            Ipv4InterfaceContainer ifc = ipv4.Assign(ndc);

            // Track connections with their addresses
            NodeConnection conn;
            conn.node1 = i;
            conn.node2 = peer;
            conn.addr1 = ifc.GetAddress(0);
            conn.addr2 = ifc.GetAddress(1);
            connections.push_back(conn);

            peerConnections[i].push_back(peer);
            peerConnections[peer].push_back(i);

            NS_LOG_INFO("Created additional connection between Node " << i << " and Node " << peer
                    << " with addresses " << conn.addr1 << " and " << conn.addr2);
        }
    }

    // Install receivers on each node and prepare peer address lists
    std::vector<Ptr<SharesReceiver>> receivers(nNodes);
    std::vector<std::vector<Ipv4Address>> peerAddresses(nNodes);
    std::vector<std::vector<uint16_t>> peerPorts(nNodes);

    for (uint32_t i = 0; i < nNodes; i++)
    {
        // Create and install the receiver application
        receivers[i] = CreateObject<SharesReceiver>();
        receivers[i]->Setup(basePort + i, i);
        nodes.Get(i)->AddApplication(receivers[i]);
        receivers[i]->SetStartTime(Seconds(0.0));
        receivers[i]->SetStopTime(Seconds(simTime));
        
        // Build the peer address list using our tracked connections
        for (uint32_t p = 0; p < peerConnections[i].size(); p++)
        {
            uint32_t peerId = peerConnections[i][p];
            
            // Find the connection between i and peerId
            for (const auto& conn : connections)
            {
                if ((conn.node1 == i && conn.node2 == peerId) ||
                    (conn.node1 == peerId && conn.node2 == i))
                {
                    // If i is node1, then get address2, else get address1
                    if (conn.node1 == i)
                    {
                        peerAddresses[i].push_back(conn.addr2);
                    }
                    else
                    {
                        peerAddresses[i].push_back(conn.addr1);
                    }
                    
                    peerPorts[i].push_back(basePort + peerId);
                    break;
                }
            }
        }
        
        NS_LOG_INFO("Node " << i << " has " << peerAddresses[i].size() 
                << " peer addresses for " << peerConnections[i].size() << " connections");
    }

    // Install sender applications on each node
    std::vector<Ptr<SharesApp>> senders(nNodes);

    for (uint32_t i = 0; i < nNodes; i++)
    {
        if (peerAddresses[i].size() > 0)
        {
            // Create and install the sender application
            senders[i] = CreateObject<SharesApp>();
            senders[i]->Setup(i, peerAddresses[i], peerPorts[i]);
            senders[i]->SetSendInterval(Seconds(sendInterval));
            senders[i]->SetRunningTime(Seconds(simTime));
            nodes.Get(i)->AddApplication(senders[i]);
            senders[i]->SetStartTime(Seconds(1.0));
            senders[i]->SetStopTime(Seconds(simTime));

            NS_LOG_INFO("Node " << i << " configured to send shares to " << peerAddresses[i].size()
                                << " peers");
        }
        else
        {
            NS_LOG_ERROR("Node " << i << " has no peer addresses configured!");
        }
    }

    // Set up position allocation for visualization
    MobilityHelper mobility;
    mobility.SetPositionAllocator(
        "ns3::GridPositionAllocator",
        "MinX",
        DoubleValue(0.0),
        "MinY",
        DoubleValue(0.0),
        "DeltaX",
        DoubleValue(nodeSpacing),
        "DeltaY",
        DoubleValue(nodeSpacing),
        "GridWidth",
        UintegerValue(static_cast<uint32_t>(std::ceil(std::sqrt(nNodes)))),
        "LayoutType",
        StringValue("RowFirst"));
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(nodes);

    // Enable packet routing
    Ipv4GlobalRoutingHelper::PopulateRoutingTables();
    
    // Create visualization
    AnimationInterface anim("mesh-p2p-shares.xml");

    // Run simulation with slightly longer time to ensure all events complete
    NS_LOG_INFO("Running simulation for " << simTime << " seconds");
    Simulator::Stop(Seconds(simTime + 1.0));  // Add 1 second buffer
    Simulator::Run();
    Simulator::Destroy();

    NS_LOG_INFO("Simulation completed successfully.");

    return 0;
}