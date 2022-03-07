#include "FractureContext.h"
#include "FractureUtil.h"

#include <Mathematics/Delaunay3.h>
#include <Mathematics/SymmetricEigensolver3x3.h>
#include <Mathematics/DistPointHyperplane.h>
#include <Mathematics/IntrLine3Plane3.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>

using namespace Deformation;

////////////////////////////////////////////////////////////////////////////////

namespace Deformation
{
	struct EdgeIntersection
	{
		// Interpreted as from mEdgeId[0] to mEdgeId[1]
		// so that a value of 0 corresponds to mEdgeId[0]
		double mParam = 0;
		glm::ivec2 mEdgeId;
	};

	// unsigned angle, our edges don't have a direction because we sort based on vertex index
	double AngleEdgePlane(const glm::dvec3& edgeDir, const glm::dvec3& planeNormal)
	{
		auto d = glm::dot(planeNormal, edgeDir);
		auto projectedEdgeDir = glm::normalize(edgeDir - d * planeNormal);
		return std::abs(glm::acos(glm::dot(projectedEdgeDir, edgeDir)));
	}

	size_t FindFaceNeighbor(
		const std::unordered_map<size_t, Tetrahedra>& tetrahedra,
		size_t mainTetId,
		const glm::ivec3& faceId
	)
	{
		for (const auto& pair : tetrahedra)
		{
			auto tetId = pair.first;
			if (tetId == mainTetId)
				continue;

			const auto& tet = pair.second;
			if (tet.ContainsFaceIndex(faceId))
				return tetId;
		}

		return -1;
	}

	std::vector<size_t> FindEdgeNeighbors(
		const std::unordered_map<size_t, Tetrahedra>& tetrahedra,
		const glm::ivec2& edgeId
	)
	{
		std::vector<size_t> neighborIds;

		for (const auto& pair : tetrahedra)
		{
			auto tetId = pair.first;
			const auto& tet = pair.second;
			if (tet.ContainsEdgeIndex(edgeId))
				neighborIds.push_back(tetId);
		}

		return neighborIds;
	}

	enum class FractureOptions
	{
		ePositiveSide,
		eNegativeSide,
		eFracture
	};

	// If this returns false we will create a new node
	// instead of snapping to the existing one
	// This should only return false if there's another tet that will be assigned ot the + side or
	// that the fracture plane will split the tet.
	bool RejectFractureDueToSpatialTolerances(
		const std::unordered_set<size_t>& neighborTets,
		const std::unordered_map<size_t, Tetrahedra>& idToTetrahedra,
		const std::vector<Vertex>& vertices,
		const glm::vec3& fracturePlanePosition,
		const glm::vec3& fracturePlaneNormal)
	{
		bool seenPositive = false;
		bool seenNegative = false;

		// So if all nodes are within this distance to the plane,
		// then we aren't going to fracture.
		// So making this value larger means we fracture less often.
		const double tolerance = 0.3;

		for (const auto& neighborId : neighborTets)
		{
			const auto& tet = idToTetrahedra.at(neighborId);

			for (const auto& nodeId : tet.mIndices)
			{
				auto d = DistPointPlane(vertices[nodeId].mPosition, fracturePlaneNormal, fracturePlanePosition);

				if (d > tolerance)
					seenPositive = true;
				if (d < -tolerance)
					seenNegative = true;

				if (seenPositive && seenNegative)
					return false;
			}
		}

		std::cout << "Returning true for the no positive node edge case." << std::endl;
		return true;
	}

	bool FractureContext::Fracture()
	{
		const auto& neighbors = GetTetrahedraNeighbors(mFractureNodeIdx);
		std::unordered_set<size_t> neighborSet(neighbors.begin(), neighbors.end());

		const double cDistanceTolerance = 0.2;
		const double cAngularSeparation = 0.34; // radians? double check what values paper uses

		// it would be nice to be able to set these tolerances to large numbers to force always snapping
		// i don't think we have this capabilitiy atm

		std::vector<EdgeIntersection> edgeIntersections;

		size_t negativeFractureNodeId = -1;

		std::array<double, 4> nodeDistToPlane;
		std::array<double, 3> edgeAngleToPlane; // really only care about the 3 edges connected to the fracture node
		auto inputStateFunctor = [&](const Tetrahedra& tet, size_t fractureNodeId)
		{
			const auto& edges = tet.GetEdges();

			for (size_t i = 0; i < 4; i++)
				nodeDistToPlane[i] = DistPointPlane(mVertices[tet.mIndices[i]].mPosition, mFracturePlaneNormal, mFractureNodePosition);

			const auto& otherNodes = tet.GetOtherVertices(fractureNodeId);

			for (size_t i = 0; i < 3; i++)
			{
				auto edgeId = GetEdgeId(fractureNodeId, otherNodes[i]);
				edgeAngleToPlane[i] = AngleEdgePlane(glm::normalize(mVertices[edgeId[1]].mPosition - mVertices[edgeId[0]].mPosition), mFracturePlaneNormal);
			}
		};

		// we won't always fracture, so the negative node is computed lazily
		auto getNegativeFractureNodeId = [&]()
		{
			if (negativeFractureNodeId == -1)
				negativeFractureNodeId = CloneVertex(mFractureNodeIdx);

			return negativeFractureNodeId;
		};

		// assign to +
		// assign to -
		// fracture (snap to edge or regular)

		bool snapToEdge = false;
		glm::ivec2 snapEdgeId;
		bool snapToFace = false;
		glm::ivec3 snapFaceId;

		auto checkTolerances = [&](const Tetrahedra& tet, size_t fractureNodeId)
		{
			{
				const float cMinVolume = 0.01f;
				if (tet.mVolume < cMinVolume)
				{
					// Experimenting with just assigning the tet to one side or the other
					const auto& otherNodes = tet.GetOtherVertices(fractureNodeId);
					snapToFace = true;
					snapFaceId = GetFaceId(fractureNodeId, fractureNodeId, otherNodes[0]);
					return;
				}
			}

			size_t numDistanceTriggered = 0;
			std::vector<size_t> closeVertexIndices;

			for (size_t i = 0; i < nodeDistToPlane.size(); i++)
			{
				if (tet.mIndices[i] == fractureNodeId)
					continue;

				if (std::abs(nodeDistToPlane[i]) <= cDistanceTolerance)
				{
					numDistanceTriggered++;
					closeVertexIndices.push_back(tet.mIndices[i]);
				}
			}

			if (numDistanceTriggered != 0)
			{
				if (numDistanceTriggered == 1)
				{
					snapToEdge = true;
					snapEdgeId = GetEdgeId(fractureNodeId, closeVertexIndices[0]);
					return;
				}
				else if (numDistanceTriggered == 2)
				{
					snapToFace = true;
					snapFaceId = GetFaceId(fractureNodeId, closeVertexIndices[0], closeVertexIndices[1]);
					return;
				}
				else
				{
					std::cout << "Small tet!" << std::endl;
					return;
				}
			}

			const auto& otherNodes = tet.GetOtherVertices(fractureNodeId);
			std::vector<size_t> snappedVertexIds;

			for (size_t i = 0; i < 3; i++)
			{
				auto angle = edgeAngleToPlane[i];
				if (angle < cAngularSeparation || angle >(glm::pi<float>() - cAngularSeparation))
					snappedVertexIds.push_back(otherNodes[i]);
			}

			if (snappedVertexIds.empty())
			{
				// this is the no snap, regular fracture case
				return;
			}
			else if (snappedVertexIds.size() == 1)
			{
				snapToEdge = true;
				snapEdgeId = GetEdgeId(mFractureNodeIdx, snappedVertexIds[0]);
				return;
			}
			else if (snappedVertexIds.size() == 2)
			{
				snapToFace = true;
				snapFaceId = GetFaceId(fractureNodeId, snappedVertexIds[0], snappedVertexIds[1]);
				return;
			}
			else
			{
				// poorly conditioned tet
				std::cout << "Very poorly conditioned tet as determined by angular threshold." << std::endl;
			}
		};

		auto onPositiveSide = [&](const Tetrahedra& tet)
		{
			glm::dvec3 center;
			for (auto nodeId : tet.mIndices)
				center += 0.25 * mVertices[nodeId].mPosition;

			return (DistPointPlane(center, mFracturePlaneNormal, mFractureNodePosition) > 0);
		};

		auto determineFractureOption = [&]()
		{
			size_t numPos = 0;
			size_t numNeg = 0;

			for (const auto& dist : nodeDistToPlane)
			{
				if (dist <= 1e-5)
					numNeg++;
				if (dist >= -1e-5)
					numPos++;
			}

			if (numPos == 4)
				return FractureOptions::ePositiveSide;
			if (numNeg == 4)
				return FractureOptions::eNegativeSide;

			return FractureOptions::eFracture;
		};

		if (RejectFractureDueToSpatialTolerances(neighborSet, mIdToTetrahedra, mVertices, mFractureNodePosition, mFracturePlaneNormal))
		{
			std::cout << "Rejecting fracture due to spatial tolerances." << std::endl;
			return false;
		}

		while (!neighborSet.empty())
		{
			// Clear shared state
			snapToEdge = false;
			snapToFace = false;
			edgeIntersections.clear();

			std::cout << "Num neighbors : " << neighborSet.size() << std::endl;

			auto fracturingTetId = *neighborSet.begin();
			const auto& tet = mIdToTetrahedra[fracturingTetId];
			inputStateFunctor(tet, mFractureNodeIdx);
			checkTolerances(tet, mFractureNodeIdx);

			if (snapToFace)
			{
				// split
				auto onPosSide = onPositiveSide(tet);
				if (onPosSide)
				{
					mNewTetrahedra.push_back(tet);
					mTetrahedraIdsToDelete.insert(fracturingTetId);
				}
				else 
				{
					auto negativeFractureNodeIdx = getNegativeFractureNodeId();

					// Assign Tet to Side
					auto copy = tet;
					copy.ReplaceVertex(mFractureNodeIdx, negativeFractureNodeIdx);
					mNewTetrahedra.push_back(copy);
					mTetrahedraIdsToDelete.insert(fracturingTetId);
				}
			}
			else if (snapToEdge)
			{
				auto edge = glm::normalize(mVertices[snapEdgeId[0]].mPosition - mVertices[snapEdgeId[1]].mPosition);
				mFracturePlaneNormal = glm::normalize(mFracturePlaneNormal - glm::dot(mFracturePlaneNormal, edge) * edge);

				inputStateFunctor(tet, mFractureNodeIdx);
			}
			//else if (snapToFace)
			//{
			//	const auto& p0 = mVertices[snapFaceId.x].mPosition;
			//	const auto& p1 = mVertices[snapFaceId.y].mPosition;
			//	const auto& p2 = mVertices[snapFaceId.z].mPosition;

			//	auto faceNormal = glm::normalize(glm::cross(p2 - p0, p1 - p0));

			//	if (glm::dot(faceNormal, mFracturePlaneNormal) < 0)
			//		faceNormal *= -1.0;

			//	mFracturePlaneNormal = faceNormal;
			//	inputStateFunctor(tet, mFractureNodeIdx);
			//}

			// Start fracturing

			// Need a bit of a restructure again,
			// where we can assign a tet to either side of the fracture plane if the tet is quite small

			// Check whether we want to fracture this tet
			auto option = determineFractureOption();
			if (option == FractureOptions::ePositiveSide)
			{
				std::cout << "Assigning to positive side" << std::endl;
				// Assign Tet to Side
				mNewTetrahedra.push_back(tet);
				mTetrahedraIdsToDelete.insert(fracturingTetId);
			}
			else if (option == FractureOptions::eNegativeSide)
			{
				std::cout << "Assigning to negative side" << std::endl;

				auto negativeFractureNodeIdx = getNegativeFractureNodeId();

				// Assign Tet to Side
				auto copy = tet;
				copy.ReplaceVertex(mFractureNodeIdx, negativeFractureNodeIdx);
				mNewTetrahedra.push_back(copy);
				mTetrahedraIdsToDelete.insert(fracturingTetId);
			}
			else
			{
				// if snap to edge, only intersect the one other edge

				if (snapToEdge)
				{
					const auto& otherNodes = tet.GetOtherVertices(mFractureNodeIdx, snapEdgeId[0] == mFractureNodeIdx ? snapEdgeId[1] : snapEdgeId[0]);
					auto edgeId = GetEdgeId(otherNodes[0], otherNodes[1]);
					double d;
					if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
						edgeIntersections.push_back({ d, edgeId });
				}
				else
				{
					// Intersect edges opposite to fracture node
					const auto& otherNodes = tet.GetOtherVertices(mFractureNodeIdx);

					{
						auto edgeId = GetEdgeId(otherNodes[0], otherNodes[1]);
						double d;
						if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
							edgeIntersections.push_back({ d, edgeId });
					}
					{
						auto edgeId = GetEdgeId(otherNodes[0], otherNodes[2]);
						double d;
						if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
							edgeIntersections.push_back({ d, edgeId });
					}
					{
						auto edgeId = GetEdgeId(otherNodes[1], otherNodes[2]);
						double d;
						if (PlaneIntersectEdge(mFractureNodePosition, mFracturePlaneNormal, mVertices[edgeId[0]].mPosition, mVertices[edgeId[1]].mPosition, d, nullptr))
							edgeIntersections.push_back({ d, edgeId });
					}
				}

				std::cout << "Num edge intersections : " << edgeIntersections.size() << std::endl;

				if (edgeIntersections.size() == 0)
				{
					// failed to assign tet to side
					throw std::exception("Failed to assign tet to side despite plane not intersecting the tet.");
				}
				else if (edgeIntersections.size() == 1)
				{
					if (!snapToEdge)
						throw std::exception("One split edge, expected edge snap but it was false.");

					auto negativeFractureNodeIdx = getNegativeFractureNodeId();

					// edge snapped fracture
					const auto& intersection = edgeIntersections.front();
					const auto& splitEdgeId = intersection.mEdgeId;
					size_t positiveSplitEdgeNode, negativeSplitEdgeNode;
					SeparateNodes(positiveSplitEdgeNode, negativeSplitEdgeNode, { (size_t)splitEdgeId[0], (size_t)splitEdgeId[1] }, mFractureNodePosition, mFracturePlaneNormal);

					size_t newEdgeNode = CreateEdgeVertex(splitEdgeId, intersection.mParam);
					size_t snappedNodeId = (snapEdgeId[0] == mFractureNodeIdx ? snapEdgeId[1] : snapEdgeId[0]);

					mTetrahedraIdsToDelete.insert(fracturingTetId);

					if (IsIsolatedEdge(splitEdgeId))
					{
						size_t positiveNewEdgeNode = newEdgeNode;
						size_t negativeNewEdgeNode = CloneVertex(positiveNewEdgeNode);

						mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, positiveNewEdgeNode, positiveSplitEdgeNode));
						mNewTetrahedra.push_back(Tetrahedra(negativeFractureNodeIdx, snappedNodeId, negativeNewEdgeNode, negativeSplitEdgeNode));
					}
					else
					{
						mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, snappedNodeId, newEdgeNode, positiveSplitEdgeNode));
						mNewTetrahedra.push_back(Tetrahedra(negativeFractureNodeIdx, snappedNodeId, newEdgeNode, negativeSplitEdgeNode));

						// case 1 - the face neighbor does not contain the snapped node, but does contain the fracture node
						// this face contains the fracture node, and both nodes from the split edge
						// for this case, tet0 contains (mFractureNodeIdx, the new edge node, the positive edge node, and the other node in the neighbor tet)
						//              , tet1 contains (negative fracture node, the new edge node, the negative edge node, and the other node in the neighbor tet)
						auto faceId0 = GetFaceId(mFractureNodeIdx, splitEdgeId[0], splitEdgeId[1]);
						auto faceNeighborId0 = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, faceId0);
						if (faceNeighborId0 != -1)
						{
							const auto& neighborTet = mIdToTetrahedra.at(faceNeighborId0);
							mTetrahedraIdsToDelete.insert(faceNeighborId0);
							mNewTetrahedra.push_back(Tetrahedra(mFractureNodeIdx, newEdgeNode, positiveSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId0[0], (size_t)faceId0[1], (size_t)faceId0[2] })));
							mNewTetrahedra.push_back(Tetrahedra(negativeFractureNodeIdx, newEdgeNode, negativeSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId0[0], (size_t)faceId0[1], (size_t)faceId0[2] })));
						}

						// case 2 - the face neighbor does not contain the fracture node, but does contain the snapped node
						// this face contains the snapped node, and both nodes from the split edge
						// for this case, tet0 contains (snapped edge node, the new edge node, the positive edge node, and the other node in the neighbor tet)
						//              , tet1 contains (snapped edge node, the new edge node, the negative edge node, and the other node in the neighbor tet)
						auto faceId1 = GetFaceId(snappedNodeId, splitEdgeId[0], splitEdgeId[1]);
						auto faceNeighborId1 = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, faceId1);
						if (faceNeighborId1 != -1)
						{
							const auto& neighborTet = mIdToTetrahedra.at(faceNeighborId1);
							mTetrahedraIdsToDelete.insert(faceNeighborId1);
							mNewTetrahedra.push_back(Tetrahedra(snappedNodeId, newEdgeNode, positiveSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId1[0], (size_t)faceId1[1], (size_t)faceId1[2] })));
							mNewTetrahedra.push_back(Tetrahedra(snappedNodeId, newEdgeNode, negativeSplitEdgeNode, neighborTet.GetOtherVertex({ (size_t)faceId1[0], (size_t)faceId1[1], (size_t)faceId1[2] })));
						}

						const auto& edgeNeighbors = FindEdgeNeighbors(mIdToTetrahedra, splitEdgeId);

						for (size_t neighborId : edgeNeighbors)
						{
							if (neighborId == fracturingTetId || neighborId == faceNeighborId0 || neighborId == faceNeighborId1)
								continue;

							const auto& edgeNeighorTet = mIdToTetrahedra[neighborId];
							const auto& remainingNodes = edgeNeighorTet.GetOtherVertices(splitEdgeId[0], splitEdgeId[1]);
							size_t posEdgeNode;
							size_t negEdgeNode;
							SeparateNodes(posEdgeNode, negEdgeNode, { (size_t)splitEdgeId[0], (size_t)splitEdgeId[1] }, mFractureNodePosition, mFracturePlaneNormal);

							mTetrahedraIdsToDelete.insert(neighborId);

							mNewTetrahedra.push_back(Tetrahedra(newEdgeNode, remainingNodes[0], remainingNodes[1], posEdgeNode));
							mNewTetrahedra.push_back(Tetrahedra(newEdgeNode, remainingNodes[0], remainingNodes[1], negEdgeNode));
						}
					}
				}
				else if (edgeIntersections.size() == 2)
				{
					const auto& otherNodes = tet.GetOtherVertices(mFractureNodeIdx);

					// regular fracture
					auto negativeFractureNodeIdx = getNegativeFractureNodeId();

					const auto& edgeId0 = edgeIntersections[0].mEdgeId;
					const auto& edgeId1 = edgeIntersections[1].mEdgeId;

					// Create two new nodes
					auto edgeVertex0 = CreateEdgeVertex(edgeId0, edgeIntersections[0].mParam);
					auto edgeVertex1 = CreateEdgeVertex(edgeId1, edgeIntersections[1].mParam);

					std::unordered_set<size_t> nonEdgeNeighbors{ fracturingTetId };

					// This remains just for validation, can remove later
					std::vector<size_t> positiveVertices;
					std::vector<size_t> negativeVertices;
					SeparateNodes(positiveVertices, negativeVertices, { otherNodes[0], otherNodes[1], otherNodes[2] }, mFractureNodePosition, mFracturePlaneNormal);

					if (!
						((positiveVertices.size() == 1 && negativeVertices.size() == 2)
							||
							(positiveVertices.size() == 2 && negativeVertices.size() == 1))
						)
					{
						throw std::exception("Expected 2 positive, 1 negative or 1 positive, 2 negative vertices during regular fracture.");
					}

					// 7 unique nodes for the regular fracture case
					// + fracture node
					// - fracture node
					// the node common to the two split edges
					// new edge node 0
					// the tetrahedral node that is on the same edge as the new edge node 0 and is not the common node
					// new edge node 1
					// the tetrahedral node ... new edge node 1 ... not common node
					//
					// if the common node is on the + side of the fracture plane, then we're in 'config 1' if not, we're in 'config 2'
					// config 1 = 1 positive node , 2 negative
					// config 2 = 2 positive nodes, 1 negative

					auto posFractureNode = mFractureNodeIdx;
					auto negFractureNode = getNegativeFractureNodeId();
					auto commonNode = GetCommonVertexFromEdges(edgeId0, edgeId1);
					auto newEdgeNode0 = edgeVertex0;
					auto edge0Node = (commonNode == edgeId0[0] ? edgeId0[1] : edgeId0[0]);
					auto newEdgeNode1 = edgeVertex1;
					auto edge1Node = (commonNode == edgeId1[0] ? edgeId1[1] : edgeId1[0]);

					// Only for isolated edges, doesn't affect neighbors
					size_t negNewEdgeNode0 = newEdgeNode0;
					size_t negNewEdgeNode1 = newEdgeNode1;

					if (IsIsolatedEdge(edgeId0))
						negNewEdgeNode0 = CloneVertex(newEdgeNode0);

					if (IsIsolatedEdge(edgeId1))
						negNewEdgeNode1 = CloneVertex(newEdgeNode1);

					auto d = DistPointPlane(mVertices[commonNode].mPosition, mFracturePlaneNormal, mFractureNodePosition);
					if (d < 0)
					{
						if (!(negativeVertices.size() == 1 && positiveVertices.size() == 2))
							throw std::exception("Expected 2 positive, 1 negative vertices during regular fracture.");

						std::swap(posFractureNode, negFractureNode);
						std::swap(newEdgeNode0, negNewEdgeNode0);
						std::swap(newEdgeNode1, negNewEdgeNode1);
					}

					// We choose the (newEdgeNode0, edge1Node) to be the diagonal edge
					mTetrahedraIdsToDelete.insert(fracturingTetId);
					mNewTetrahedra.push_back(Tetrahedra(posFractureNode, commonNode, newEdgeNode0, newEdgeNode1));
					mNewTetrahedra.push_back(Tetrahedra(negFractureNode, negNewEdgeNode0, edge0Node, edge1Node));
					mNewTetrahedra.push_back(Tetrahedra(negFractureNode, negNewEdgeNode0, negNewEdgeNode1, edge1Node));

					// the first face neighbor, containing the fracture node, and both nodes from a single split edge
					{
						auto fracturedFaceId = GetFaceId(mFractureNodeIdx, edgeId0[0], edgeId0[1]);
						auto faceNeighborId = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, fracturedFaceId);
						if (faceNeighborId != -1)
						{
							const auto& faceNeighborTet = mIdToTetrahedra[faceNeighborId];
							mNewTetrahedra.push_back(Tetrahedra(posFractureNode, newEdgeNode0, commonNode, faceNeighborTet.GetOtherVertex({ (size_t)mFractureNodeIdx, (size_t)edgeId0[0], (size_t)edgeId0[1] })));
							mNewTetrahedra.push_back(Tetrahedra(negFractureNode, newEdgeNode0, edge0Node, faceNeighborTet.GetOtherVertex({ (size_t)mFractureNodeIdx, (size_t)edgeId0[0], (size_t)edgeId0[1] })));
							mTetrahedraIdsToDelete.insert(faceNeighborId);
							nonEdgeNeighbors.insert(faceNeighborId);
						}
					}
					// the second face neighbor, containing the fracture node, and both nodes from the other split edge
					{
						auto fracturedFaceId = GetFaceId(mFractureNodeIdx, edgeId1[0], edgeId1[1]);
						auto faceNeighborId = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, fracturedFaceId);
						if (faceNeighborId != -1)
						{
							const auto& faceNeighborTet = mIdToTetrahedra[faceNeighborId];
							mNewTetrahedra.push_back(Tetrahedra(posFractureNode, newEdgeNode1, commonNode, faceNeighborTet.GetOtherVertex({ (size_t)mFractureNodeIdx, (size_t)edgeId1[0], (size_t)edgeId1[1] })));
							mNewTetrahedra.push_back(Tetrahedra(negFractureNode, newEdgeNode1, edge1Node, faceNeighborTet.GetOtherVertex({ (size_t)mFractureNodeIdx, (size_t)edgeId1[0], (size_t)edgeId1[1] })));
							mTetrahedraIdsToDelete.insert(faceNeighborId);
							nonEdgeNeighbors.insert(faceNeighborId);
						}
					}
					// the third face neighbor, containing the nodes that are not the fracture node
					{
						auto fracturedFaceId = GetFaceId(commonNode, edge0Node, edge1Node);
						auto faceNeighborId = FindFaceNeighbor(mIdToTetrahedra, fracturingTetId, fracturedFaceId);
						if (faceNeighborId != -1)
						{
							const auto& faceNeighborTet = mIdToTetrahedra[faceNeighborId];
							auto otherNode = faceNeighborTet.GetOtherVertex({ (size_t)commonNode, (size_t)edge0Node, (size_t)edge1Node });
							mNewTetrahedra.push_back(Tetrahedra(commonNode, newEdgeNode0, newEdgeNode1, otherNode));
							// newEdgeNode0 and edge1Node form the diagonal
							mNewTetrahedra.push_back(Tetrahedra(edge0Node, newEdgeNode0, edge1Node, otherNode));
							mNewTetrahedra.push_back(Tetrahedra(newEdgeNode1, newEdgeNode0, edge1Node, otherNode));
							mTetrahedraIdsToDelete.insert(faceNeighborId);
							nonEdgeNeighbors.insert(faceNeighborId);
						}
					}

					// TODO think about whether we can remove the separate nodes below in favor of using the previously defined nodes above
					{
						const auto& edgeNeighbors = FindEdgeNeighbors(mIdToTetrahedra, edgeId0);

						for (size_t neighborId : edgeNeighbors)
						{
							if (nonEdgeNeighbors.find(neighborId) != nonEdgeNeighbors.end())
								continue;

							const auto& edgeNeighorTet = mIdToTetrahedra[neighborId];
							const auto& remainingNodes = edgeNeighorTet.GetOtherVertices(edgeId0[0], edgeId0[1]);
							size_t posEdgeNode;
							size_t negEdgeNode;
							SeparateNodes(posEdgeNode, negEdgeNode, { (size_t)edgeId0[0], (size_t)edgeId0[1] }, mFractureNodePosition, mFracturePlaneNormal);

							mTetrahedraIdsToDelete.insert(neighborId);

							mNewTetrahedra.push_back(Tetrahedra(edgeVertex0, remainingNodes[0], remainingNodes[1], posEdgeNode));
							mNewTetrahedra.push_back(Tetrahedra(edgeVertex0, remainingNodes[0], remainingNodes[1], negEdgeNode));
						}
					}
					{
						const auto& edgeNeighbors = FindEdgeNeighbors(mIdToTetrahedra, edgeId1);

						for (size_t neighborId : edgeNeighbors)
						{
							if (nonEdgeNeighbors.find(neighborId) != nonEdgeNeighbors.end())
								continue;

							const auto& edgeNeighorTet = mIdToTetrahedra[neighborId];
							const auto& remainingNodes = edgeNeighorTet.GetOtherVertices(edgeId1[0], edgeId1[1]);
							size_t posEdgeNode;
							size_t negEdgeNode;
							SeparateNodes(posEdgeNode, negEdgeNode, { (size_t)edgeId1[0], (size_t)edgeId1[1] }, mFractureNodePosition, mFracturePlaneNormal);

							mTetrahedraIdsToDelete.insert(neighborId);

							mNewTetrahedra.push_back(Tetrahedra(edgeVertex1, remainingNodes[0], remainingNodes[1], posEdgeNode));
							mNewTetrahedra.push_back(Tetrahedra(edgeVertex1, remainingNodes[0], remainingNodes[1], negEdgeNode));
						}
					}
				}
				else
				{
					// tetrahedra is small
				}
			}

			for (auto tetId : mTetrahedraIdsToDelete)
			{
				mIdToTetrahedra.erase(tetId);
				neighborSet.erase(tetId);
			}

			for (auto& tet : mNewTetrahedra)
				mIdToTetrahedra[mTetIdCounter++] = tet;

			// verify that all nodes are referenced by at least one tet
			{
				std::unordered_set<size_t> vertexIds;

				for (const auto& pair : mIdToTetrahedra)
				{
					for (const auto& vertexId : pair.second.mIndices)
						vertexIds.insert(vertexId);
				}

				for (size_t id = 0; id < mVertices.size(); id++)
				{
					if (vertexIds.find(id) == vertexIds.end())
					{
						std::cout << "Vertex id : " << id << std::endl;
						throw std::exception("Isolated vertex detected.");
					}
				}
			}

			mTetrahedraIdsToDelete.clear();
			mNewTetrahedra.clear();
		}

		return true;
	}

	////////////////////////////////////////////////////////////////////////////////

	size_t FractureContext::CloneVertex(size_t vertexId)
	{
		Vertex vertex;
		vertex.mPosition = mVertices[vertexId].mPosition;
		vertex.mMaterialCoordinates = mVertices[vertexId].mMaterialCoordinates;
		vertex.mVelocity = mVertices[vertexId].mVelocity;

		mVertices.push_back(vertex);

		return (mVertices.size() - 1);
	}

	////////////////////////////////////////////////////////////////////////////////

	size_t FractureContext::CreateEdgeVertex(const glm::ivec2& edgeId, double parametricDistance)
	{
		Vertex vertex;
		const auto& p0 = mVertices[edgeId.x].mPosition;
		const auto& p1 = mVertices[edgeId.y].mPosition;

		vertex.mPosition = p0 + (p1 - p0) * parametricDistance;

		const auto& m0 = mVertices[edgeId.x].mMaterialCoordinates;
		const auto& m1 = mVertices[edgeId.y].mMaterialCoordinates;
		vertex.mMaterialCoordinates = m0 + (m1 - m0) * parametricDistance;

		const auto& v0 = mVertices[edgeId.x].mVelocity;
		const auto& v1 = mVertices[edgeId.y].mVelocity;
		// Not sure this is the right way to do this but it seems rather reasonable
		vertex.mVelocity = v0 + (v1 - v0) * parametricDistance;

		mVertices.push_back(vertex);

		return (mVertices.size() - 1);
	}

	////////////////////////////////////////////////////////////////////////////////

	void FractureContext::SeparateNodes(size_t& outPositiveNode, size_t& outNegativeNode, const std::array<size_t, 2>& nodes, const glm::dvec3& planePos, const glm::dvec3& planeNormal) const
	{
		auto d0 = DistPointPlane(mVertices[nodes[0]].mPosition, planeNormal, planePos);
		auto d1 = DistPointPlane(mVertices[nodes[1]].mPosition, planeNormal, planePos);

		size_t numPos = 0;
		size_t numNeg = 0;

		if (d0 < 0)
			numNeg++;
		if (d1 < 0)
			numNeg++;
		if (d0 > 0)
			numPos++;
		if (d1 > 0)
			numPos++;

		if (!(numPos == 1 && numNeg == 1))
		{
			std::cout << "Num Pos : " << numPos << " Num Neg : " << numNeg << std::endl;
			std::cout << "d0 : " << d0 << " d1 : " << d1 << std::endl;
			throw std::exception("Expected one positive and one negative vertex in isolated edge fracture.");
		}

		outPositiveNode = (d0 > 0 ? nodes[0] : nodes[1]);
		outNegativeNode = (d0 < 0 ? nodes[0] : nodes[1]);
	}

	////////////////////////////////////////////////////////////////////////////////

	void FractureContext::SeparateNodes(std::vector<size_t>& outPositiveNodes, std::vector<size_t>& outNegativeNodes, const std::vector<size_t>& nodes, const glm::dvec3& planePos, const glm::dvec3& planeNormal) const
	{
		for (size_t idx : nodes)
		{
			auto d = DistPointPlane(mVertices[idx].mPosition, planeNormal, planePos);

			if (d > 0)
				outPositiveNodes.push_back(idx);
			else if (d < 0)
				outNegativeNodes.push_back(idx);
		}
	}

	////////////////////////////////////////////////////////////////////////////////

	bool FractureContext::PlaneIntersectEdge(const glm::dvec3& planePos, const glm::dvec3& planeNormal, const glm::dvec3& edgePos0, const glm::dvec3& edgePos1, double& d, glm::dvec3* intersectionPos) const
	{
		gte::Vector3<double> normal
		{
			planeNormal.x,
			planeNormal.y,
			planeNormal.z
		};
		gte::Vector3<double> planeOrigin
		{
		   planePos.x,
		   planePos.y,
		   planePos.z
		};
		gte::Plane3<double> plane(normal, planeOrigin);

		auto dir = glm::normalize(edgePos1 - edgePos0);

		gte::Vector3<double> lineDir
		{
			dir.x, dir.y, dir.z
		};

		gte::Vector3<double> lineOrigin
		{
			edgePos0.x, edgePos0.y, edgePos0.z
		};

		gte::Line3<double> line(lineOrigin, lineDir);

		// Get signed distance of all vertices to plane
		gte::FIQuery<double, gte::Line3<double>, gte::Plane3<double>> findIntersectionQuery;
		auto results = findIntersectionQuery(line, plane);

		if (!results.intersect)
			return false;

		if (isnan(results.parameter))
			throw std::exception("Nan detected in plane-edge intersection.");

		std::cout << "Result parameter : " << results.parameter << std::endl;

		double param = results.parameter / (glm::distance(edgePos0, edgePos1));

		std::cout << "Scaled Result parameter : " << param << std::endl;

		if (param <= 0 || param >= 1.0)
			return false;

		d = param;

		if (intersectionPos)
		{
			intersectionPos->x = results.point[0];
			intersectionPos->y = results.point[1];
			intersectionPos->z = results.point[2];
		}

		return true;
	}

	////////////////////////////////////////////////////////////////////////////////

	bool FractureContext::IsIsolatedEdge(const glm::ivec2& edgeId) const
	{
		size_t counter = 0;

		for (const auto& pair : mIdToTetrahedra)
		{
			const auto& tet = pair.second;
			if (tet.ContainsEdgeIndex(edgeId))
				counter++;
		}

		if (counter == 0)
			throw std::exception("Failed to find isolated edge id in any tetrahedra.");

		return (counter == 1);
	}

	////////////////////////////////////////////////////////////////////////////////

	std::vector<size_t> FractureContext::GetTetrahedraNeighbors(size_t nodeIdx) const
	{
		std::vector<size_t> neighborIds;

		for (const auto& pair : mIdToTetrahedra)
		{
			size_t tetId = pair.first;
			const auto& tetrahedra = pair.second;
			if (tetrahedra.ContainsVertexIndex(nodeIdx))
				neighborIds.push_back(tetId);
		}

		return neighborIds;
	}

	////////////////////////////////////////////////////////////////////////////////

	size_t FractureContext::GetNonFractureNode(const glm::ivec2& edge) const
	{
		if (!(edge[0] == mFractureNodeIdx || edge[1] == mFractureNodeIdx))
			throw std::exception("Expected fracture node to be present in GetNonFractureNode.");

		return (edge[0] == mFractureNodeIdx ? edge[1] : edge[0]);
	}
}

////////////////////////////////////////////////////////////////////////////////