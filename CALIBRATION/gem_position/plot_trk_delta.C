#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <map>

const int MAX_DATA  = 10000;
const int maxTracks = 100;
const int nPlanes   = 5;

const double d0_cut               = 20.0;
const double delta_trans_cut[2]   = {-100, 100.0};
const double delta_long_cut[2]    = {-50.0, 50.0};
const double plane_theta[nPlanes] = {150.0, 150.0, 127.0, 127.0, 104.0}; // Angle in degrees
const double plane_r[nPlanes]     = {615.0, 655.6, 523.0, 563.6, 615.0}; // Radius of the second point

// Should ultimately get these from the param file
struct GEMParams {
  double theta[2]; // GEM angles in degrees
  double phi[2];   // GEM phi angles in degrees
  double r[2];     // GEM radii
  double dx[2];    // GEM dx offsets
  double dy[2];    // GEM dy offsets
};

// Map to store GEMParams for each run number
std::map<int, GEMParams> runToGEMParams = {
    {0, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 0.0}, {0.0, 0.0}}},
    {300000, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 0.0}, {0.0, 0.0}}},
    {300001, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 0.5}, {0.0, 0.0}}},
    {300002, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.0}, {0.0, 0.0}}},
    {300003, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.3}, {0.0, 0.0}}},
    {300004, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.4}, {0.0, 0.0}}},
    {300005, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {0.0, 0.0}}},
    {300006, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.6}, {0.0, 0.0}}},
    {300007, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.7}, {0.0, 0.0}}},
    {300008, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.8}, {0.0, 0.0}}},
    {300009, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.9}, {0.0, 0.0}}},
    {300010, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 2.0}, {0.0, 0.0}}},
    {300011, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.000}, {0.0, 1.5}, {0.0, 0.0}}},
    {300012, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 96.000}, {0.0, 1.5}, {0.0, 0.0}}},

    {300020, {{127.500, 127.500}, {0.0, 0.0}, {72.5324, 90.0596}, {3.64493, 4.50118}, {0.0, 0.0}}},
    {300021, {{126.729, 126.729}, {0.0, 0.0}, {84.4906, 104.183}, {-0.909749, -1.03288}, {0.0, 0.5}}},
    {300022, {{127.000, 127.000}, {0.0, 0.0}, {77.571, 95.571}, {0.0, 1.5}, {1.5, 1.5}}},
    {300023, {{126.737, 126.737}, {0.0, 0.0}, {78.8751, 97.0595}, {-0.0971787, -0.220268}, {1.5, 1.5}}},
    {300024, {{124.410, 124.410}, {0.0, 0.0}, {67.9565, 86.4893}, {-3.97332, -5.7903}, {1.5, 1.5}}},
    {300025, {{125.425, 125.425}, {0.0, 0.0}, {85.8521, 106.21}, {-4.26412, -5.16062}, {1.5, 1.5}}},
    {300026, {{125.425, 125.425}, {0.0, 0.0}, {87.0633, 105.807}, {-7.51136, -8.1162}, {1.5, 1.5}}},
    {300027, {{127.941, 127.941}, {0.0, 0.0}, {83.714, 102.727}, {-2.81105, -3.50819}, {1.5, 1.5}}},
    {300028, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {1.5, 1.5}, {1.5, 1.5}}},
    {300029, {{127.0, 127.0}, {0.0, 0.0}, {77.571, 95.571}, {1.5, 1.5}, {1.5, 1.5}}}};

double gem_theta = 127.0;            // Angle in degrees
double gem_phi   = 0.0;              // Angle in degrees
double gem_r[2]  = {77.571, 95.571}; // Radius of the GEM's
double gem_dx[2] = {0.0, 0.0};       // GEM dx offsets
double gem_dy[2] = {0.0, 0.0};       // GEM dy offsets

const int nFixedz = 3;
// const double target_z[nFixedz] = {0.0}; // Fixed z positions for the planes
const double target_z[nFixedz] = {-10.0, 0.0, 10.0}; // Fixed z positions for the planes
// const double z_range[nFixedz + 1] = {-15.0, -5.0, 5.0, 15.0}; // Tolerance range for the fixed z positions (low,
// high)
const double z_range[nFixedz + 1] = {-20.0, -5, 5.0, 20.0}; // Tolerance range for the fixed z positions (low,
// high)
// const double z_range[nFixedz + 1] = {-20.0, 20.0}; // Tolerance range for the fixed z positions (low, high)

const bool use_projz = true; // Fix z position to GEM projz

const double dx_min = -30.0;
const double dx_max = 30.0;
const int dx_NBINS  = 60;
const double dy_min = -30.0;
const double dy_max = 30.0;
const int dy_NBINS  = 60;

TVector3 LinePlaneIntersection(const TVector3 &p1, const TVector3 &p2, const TVector3 &p3, const TVector3 &l1,
                               const TVector3 &l2) {
  // Define the plane normal
  TVector3 planeNormal = (p2 - p1).Cross(p3 - p1);
  planeNormal          = planeNormal.Unit();

  // Line direction
  TVector3 lineDir = l2 - l1;

  // Check if the line is parallel to the plane
  double denom = planeNormal.Dot(lineDir);
  if (fabs(denom) < 1e-6) {
    throw std::runtime_error("The line is parallel to the plane and does not intersect.");
  }

  // Calculate the intersection point
  double t              = planeNormal.Dot(p1 - l1) / denom;
  TVector3 intersection = l1 + t * lineDir;

  return intersection;
}

TVector3 GetHodoHitPosition(const int paddle, const int plane) {

  double dTrans = (5 - paddle) * 22.0;
  // Define the plane using three points
  TVector3 p_hit(plane_r[plane] * cos((plane_theta[plane] - 90) * TMath::DegToRad()), 0,
                 -plane_r[plane] * sin((plane_theta[plane] - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-dTrans * cos((180 - plane_theta[plane]) * TMath::DegToRad()), 0,
                           -dTrans * sin((180 - plane_theta[plane]) * TMath::DegToRad()));

  return p_hit;
}

void plot_trk_delta() {
  int run_number = 400002;
  // Get the GEM parameters for the run number
  auto it = runToGEMParams.find(run_number);
  if (it != runToGEMParams.end()) {
    gem_theta = it->second.theta[0]; // Use the first GEM's theta
    gem_r[0]  = it->second.r[0];     // Use the first GEM's radius
    gem_dx[0] = it->second.dx[0];    // Use the first GEM's dx
    gem_dy[0] = it->second.dy[0];    // Use the first GEM's dy
    gem_r[1]  = it->second.r[1];     // Use the second GEM's radius
    gem_dx[1] = it->second.dx[1];    // Use the second GEM's dx
    gem_dy[1] = it->second.dy[1];    // Use the second GEM's dy
  }

  // Define the plane using three points
  TVector3 p1[2], p2[2], p3[2];
  for (int i = 0; i < 2; ++i) {
    double theta = gem_theta;
    p1[i]        = TVector3(gem_r[i] * cos((theta - 90) * TMath::DegToRad()), 0,
                            -gem_r[i] * sin((theta - 90) * TMath::DegToRad()));
    p2[i]        = p1[i] + TVector3(-gem_r[i] * cos((180 - theta) * TMath::DegToRad()), 0,
                                    -gem_r[i] * sin((180 - theta) * TMath::DegToRad()));
    p3[i]        = p1[i] + TVector3(0, 10.0, 0); // Arbitrary y value of 10.0
  }

  ///////////////////////////////////////////////////////////////////
  TString fileName = "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/"
                     //  "LAD_COIN_22282_-1_inverted.root";
                     //  "LAD_COIN_22282_-1_500trks_good_timing.root";
                     //  "LAD_COIN_22282_" + std::to_string(run_number) + ".root";
                     "LAD_COIN_23108_4_4_-1.root";
  //  "LAD_COIN_22073_-1.root";
  TString outputFileName = "gem_pos_fixing_23108_4_4_-1.root";
  // Open the ROOT file
  TFile *file = TFile::Open(fileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open the ROOT file!" << std::endl;
    return;
  }

  // Get the TTree
  TTree *T = dynamic_cast<TTree *>(file->Get("T"));
  if (!T) {
    std::cerr << "Error: Cannot find the TTree named 'T'!" << std::endl;
    file->Close();
    return;
  }
  // Define arrays to hold the data
  Double_t trk_d0[MAX_DATA];
  Double_t trk_d0_good[MAX_DATA];
  Double_t trk_projz[MAX_DATA];
  Double_t trk_x[2][MAX_DATA], trk_y[2][MAX_DATA], trk_z[2][MAX_DATA];
  Double_t trk_x_local[2][MAX_DATA], trk_y_local[2][MAX_DATA];
  Double_t kin_trackID_0[MAX_DATA], kin_trackID_1[MAX_DATA];
  Double_t kin_plane_0[MAX_DATA], kin_plane_1[MAX_DATA];
  Double_t kin_paddle_0[MAX_DATA], kin_paddle_1[MAX_DATA];
  Double_t kin_hittime_0[MAX_DATA], kin_hittime_1[MAX_DATA];
  Double_t kin_hittheta_0[MAX_DATA], kin_hittheta_1[MAX_DATA];
  Double_t kin_hitphi_0[MAX_DATA], kin_hitphi_1[MAX_DATA];
  Double_t kin_hitedep_0[MAX_DATA], kin_hitedep_1[MAX_DATA];
  Double_t kin_dTrkHoriz_0[MAX_DATA], kin_dTrkHoriz_1[MAX_DATA];
  Double_t kin_dTrkVert_0[MAX_DATA], kin_dTrkVert_1[MAX_DATA];
  Int_t nTracks, nGoodHits;
  Double_t vertex_x, vertex_y, vertex_z;

  T->SetBranchAddress("Ndata.H.gem.trk.d0", &nTracks);
  T->SetBranchAddress("Ndata.H.ladkin.goodhit_trackid_0", &nGoodHits);
  T->SetBranchAddress("H.gem.trk.d0", &trk_d0);
  T->SetBranchAddress("H.gem.trk.d0_good", &trk_d0_good);
  T->SetBranchAddress("H.gem.trk.projz", &trk_projz);
  T->SetBranchAddress("H.gem.trk.x1", &trk_x[0]);
  T->SetBranchAddress("H.gem.trk.y1", &trk_y[0]);
  T->SetBranchAddress("H.gem.trk.z1", &trk_z[0]);
  T->SetBranchAddress("H.gem.trk.x2", &trk_x[1]);
  T->SetBranchAddress("H.gem.trk.y2", &trk_y[1]);
  T->SetBranchAddress("H.gem.trk.z2", &trk_z[1]);
  T->SetBranchAddress("H.gem.trk.x1_local", &trk_x_local[0]);
  T->SetBranchAddress("H.gem.trk.y1_local", &trk_y_local[0]);
  T->SetBranchAddress("H.gem.trk.x2_local", &trk_x_local[1]);
  T->SetBranchAddress("H.gem.trk.y2_local", &trk_y_local[1]);
  T->SetBranchAddress("H.ladkin.goodhit_trackid_0", &kin_trackID_0);
  T->SetBranchAddress("H.ladkin.goodhit_trackid_1", &kin_trackID_1);
  T->SetBranchAddress("H.ladkin.goodhit_plane_0", &kin_plane_0);
  T->SetBranchAddress("H.ladkin.goodhit_plane_1", &kin_plane_1);
  T->SetBranchAddress("H.ladkin.goodhit_paddle_0", &kin_paddle_0);
  T->SetBranchAddress("H.ladkin.goodhit_paddle_1", &kin_paddle_1);
  T->SetBranchAddress("H.ladkin.goodhit_hittime_0", &kin_hittime_0);
  T->SetBranchAddress("H.ladkin.goodhit_hittime_1", &kin_hittime_1);
  T->SetBranchAddress("H.ladkin.goodhit_hittheta_0", &kin_hittheta_0);
  T->SetBranchAddress("H.ladkin.goodhit_hittheta_1", &kin_hittheta_1);
  T->SetBranchAddress("H.ladkin.goodhit_hitphi_0", &kin_hitphi_0);
  T->SetBranchAddress("H.ladkin.goodhit_hitphi_1", &kin_hitphi_1);
  T->SetBranchAddress("H.ladkin.goodhit_hitedep_0", &kin_hitedep_0);
  T->SetBranchAddress("H.ladkin.goodhit_hitedep_1", &kin_hitedep_1);
  T->SetBranchAddress("H.ladkin.goodhit_dTrkHoriz_0", &kin_dTrkHoriz_0);
  T->SetBranchAddress("H.ladkin.goodhit_dTrkHoriz_1", &kin_dTrkHoriz_1);
  T->SetBranchAddress("H.ladkin.goodhit_dTrkVert_0", &kin_dTrkVert_0);
  T->SetBranchAddress("H.ladkin.goodhit_dTrkVert_1", &kin_dTrkVert_1);
  T->SetBranchAddress("H.react.x", &vertex_x);
  T->SetBranchAddress("H.react.y", &vertex_y);
  T->SetBranchAddress("H.react.z", &vertex_z);

  // Create a histogram to store the data
  TH2D *h_dx_gemX[2];
  TH2D *h_dx_dy_gem[2];
  TH2D *h_dx_dy_gem_good_vertex[2];
  TH2D *h_dx_dy_gem_bad_vertex[2];
  TH2D *h_gem_xz[2];
  TH2D *h_intersection_xz[2];
  TH2D *h_gem_xz_good_vertex[2];
  TH2D *h_gem_xz_bad_vertex[2];
  TH2D *h_intersection_xz_good_vertex[2];
  TH2D *h_intersection_xz_bad_vertex[2];

  for (int i = 0; i < 2; ++i) {
    h_dx_gemX[i] = new TH2D(Form("h_dx_gemX_gem%d", i), Form("dx vs GEM X for GEM %d;dx (cm);GEM X (cm)", i), dx_NBINS,
                            dx_min, dx_max, 120, -60, 60);

    h_dx_dy_gem[i] = new TH2D(Form("h_dx_dy_gem%d", i), Form("dx vs dy distribution for GEM %d;dx (cm);dy (cm)", i),
                              dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);

    h_dx_dy_gem_good_vertex[i] = new TH2D(Form("h_dx_dy_gem%d_good_vertex", i),
                                          Form("dx vs dy distribution for GEM %d (Good Vertex);dx (cm);dy (cm)", i),
                                          dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);

    h_dx_dy_gem_bad_vertex[i] = new TH2D(Form("h_dx_dy_gem%d_bad_vertex", i),
                                         Form("dx vs dy distribution for GEM %d (Bad Vertex);dx (cm);dy (cm)", i),
                                         dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);

    h_gem_xz[i] = new TH2D(Form("h_gem%d_xz", i), Form("GEM %d x vs z;x (cm);z (cm)", i), 120, 0, 120, 200, -160, 40);
    h_gem_xz_good_vertex[i] =
        new TH2D(Form("h_gem%d_xz_good_vertex", i), Form("GEM %d x vs z (Good Vertex);x (cm);z (cm)", i), 120, 0, 120,
                 200, -160, 40);
    h_gem_xz_bad_vertex[i] = new TH2D(Form("h_gem%d_xz_bad_vertex", i),
                                      Form("GEM %d x vs z (Bad Vertex);x (cm);z (cm)", i), 120, 0, 120, 200, -160, 40);
    h_intersection_xz[i]   = new TH2D(Form("h_intersection%d_xz", i), Form("Intersection %d x vs z;x (cm);z (cm)", i),
                                      120, 0, 120, 200, -160, 40);
    h_intersection_xz_good_vertex[i] =
        new TH2D(Form("h_intersection%d_xz_good_vertex", i),
                 Form("Intersection %d x vs z (Good Vertex);x (cm);z (cm)", i), 120, 0, 120, 200, -160, 40);
    h_intersection_xz_bad_vertex[i] =
        new TH2D(Form("h_intersection%d_xz_bad_vertex", i),
                 Form("Intersection %d x vs z (Bad Vertex);x (cm);z (cm)", i), 120, 0, 120, 200, -160, 40);
  }
  TH2D *h_gem1dx_vs_gem2dx = new TH2D("h_gem1dx_vs_gem2dx", "GEM1 dx vs GEM2 dx;GEM1 dx (cm);GEM2 dx (cm)", dx_NBINS,
                                      dx_min, dx_max, dx_NBINS, dx_min, dx_max);
  TH2D *h_gem1dx_vs_gem2dx_good =
      new TH2D("h_gem1dx_vs_gem2dx_good", "GEM1 dx vs GEM2 dx (Good Vertex);GEM1 dx (cm);GEM2 dx (cm)", dx_NBINS,
               dx_min, dx_max, dx_NBINS, dx_min, dx_max);
  TH2D *h_gem1dx_vs_gem2dx_bad =
      new TH2D("h_gem1dx_vs_gem2dx_bad", "GEM1 dx vs GEM2 dx (Bad Vertex);GEM1 dx (cm);GEM2 dx (cm)", dx_NBINS, dx_min,
               dx_max, dx_NBINS, dx_min, dx_max);
  TH2D *h_gem1dx_vs_plane = new TH2D("h_gem1dx_vs_plane", "GEM1 dx vs Plane;Plane;GEM1 dx (cm)", nPlanes, 0, nPlanes,
                                     dx_NBINS, dx_min, dx_max);
  TH2D *h_gem1dx_vs_paddle[nPlanes];
  for (int planeIndex = 0; planeIndex < nPlanes; ++planeIndex) {
    h_gem1dx_vs_paddle[planeIndex] = new TH2D(Form("h_gem1dx_vs_paddle_plane%d", planeIndex),
                                              Form("GEM1 dx vs Paddle for Plane %d;Paddle;GEM1 dx (cm)", planeIndex),
                                              11, -0.5, 10.5, dx_NBINS, dx_min, dx_max);
  }
  TH2D *h_hodo_hit_xz = new TH2D("h_hodo_hit_xz", "Hodoscope Hit x vs z;x (cm);z (cm)", 600, 0, 600, 600, -700, 0);
  // Create histograms for target x vs z
  TH2D *h_target_xz = new TH2D("h_target_xz", "Target x vs z;x (cm);z (cm)", 40, -4, 4, 100, -20, 20);

  TH2D *h_gem1dx_vs_targetz = new TH2D("h_gem1dx_vs_targetz", "GEM1 dx vs Target z;GEM1 dx (cm);Target z (cm)",
                                       dx_NBINS, dx_min, dx_max, 100, -20, 20);
  TH1D *h_trk_projz         = new TH1D("h_trk_projz", "Track projz;projz (cm);Counts", 100, -20, 20);
  // Create a 2D histogram for projz vs vertex z
  TH2D *h_projz_vs_vertexz =
      new TH2D("h_projz_vs_vertexz", "projz vs vertex z;projz (cm);vertex z (cm)", 100, -20, 20, 100, -20, 20);
  ///////////////////////////////////////////////////////////////////

  // Loop through the tree entries
  Long64_t nEntries = T->GetEntries();
  // nEntries          = 10000; // For testing purposes, limit to 1000 entries
  for (Long64_t i = 0; i < nEntries; ++i) {
    T->GetEntry(i);
    // Skip processing if there are too many tracks
    if (nTracks > maxTracks) {
      continue;
    }

    // Create a vector to label tracks as good or not
    std::vector<bool> goodTrack(nTracks, false);
    // Create vectors to store hodoscope paddle, plane, and y position
    std::vector<int> hodoPaddle(nTracks, -1);
    std::vector<int> hodoPlane(nTracks, -1);
    std::vector<double> hodoYPos(nTracks, 0.0);
    // Loop through all hodo hits
    for (int j = 0; j < nGoodHits; ++j) {
      // Get the track ID associated with the hodo hit
      int trackID = static_cast<int>(kin_trackID_0[j]);

      // cout << "Processing entry " << i << " with " << nTracks << " tracks and " << nGoodHits << " good hodo hits."
      //      << endl;

      // Ensure the track ID is within bounds
      if (trackID < 0 || trackID >= nTracks) {
        continue;
      }

      // Apply the cuts
      if (fabs(trk_d0[trackID]) < d0_cut &&                                                     // d0 cut
          kin_dTrkHoriz_0[j] > delta_trans_cut[0] && kin_dTrkHoriz_0[j] < delta_trans_cut[1] && // delta_pos_trans cut
          kin_dTrkVert_0[j] > delta_long_cut[0] && kin_dTrkVert_0[j] < delta_long_cut[1]) {     // delta_pos_long cut
        if (goodTrack[trackID]) {
          std::cout << "Warning: Track ID " << trackID << " is already marked as good!" << std::endl;
        }
        goodTrack[trackID]  = true;
        hodoPaddle[trackID] = static_cast<int>(kin_paddle_0[j]);
        hodoPlane[trackID]  = static_cast<int>(kin_plane_0[j]);
        // TODO: Need to get Y position somehow
      }
    }

    // Loop through all tracks
    for (int j = 0; j < nTracks; ++j) {
      // Fill the histograms based on the goodTrack vector
      if (!goodTrack[j]) {
        continue;
      }
      // Fill the histogram for projz
      h_trk_projz->Fill(trk_projz[j]);
      h_projz_vs_vertexz->Fill(trk_projz[j], vertex_z);

      TVector3 p_hit = GetHodoHitPosition(hodoPaddle[j], hodoPlane[j]);
      TVector3 p_react(vertex_x, vertex_y, vertex_z);
      bool good_vertez = true;
      if (nFixedz <= 0) {
        if (vertex_x > 1000 || vertex_y > 1000 || vertex_z > 1000) {
          good_vertez = false;
        }
        if (vertex_x > 1000)
          p_react.SetX(0);
        if (vertex_y > 1000)
          p_react.SetY(0);
        if (vertex_z > 1000)
          p_react.SetZ(0);
      } else {
        good_vertez = false;
        for (int k = 0; k < nFixedz; ++k) {
          if (use_projz) {
            if (trk_projz[j] > z_range[k] && trk_projz[j] < z_range[k + 1]) {
              good_vertez = true;
              p_react.SetZ(target_z[k]);
              break;
            }
          } else {
            if (vertex_z > z_range[k] && vertex_z < z_range[k + 1]) {
              good_vertez = true;
              p_react.SetZ(target_z[k]);
              break;
            }
          }
        }
      }

      // p_react.SetZ(-p_react.Z());
      // Calculate the intersection point

      // Fill the histograms with the calculated values
      double xz_separation[2];
      for (int gemIndex = 0; gemIndex < 2; ++gemIndex) {
        TVector3 interaction = LinePlaneIntersection(p1[gemIndex], p2[gemIndex], p3[gemIndex], p_hit, p_react);
        TVector3 p_gem_hit(trk_x[gemIndex][j], trk_y[gemIndex][j], trk_z[gemIndex][j]);
        // Calculate the x-z separation and y separation
        xz_separation[gemIndex] =
            sqrt(pow(p_gem_hit.X() - interaction.X(), 2) + pow(p_gem_hit.Z() - interaction.Z(), 2));
        double x_diff = p_gem_hit.X() - interaction.X();
        xz_separation[gemIndex] *= (x_diff > 0 ? 1 : -1);
        double y_separation = p_gem_hit.Y() - interaction.Y();

        h_dx_dy_gem[gemIndex]->Fill(xz_separation[gemIndex], y_separation);
        h_gem_xz[gemIndex]->Fill(p_gem_hit.X(), p_gem_hit.Z());

        if (good_vertez) {
          h_dx_dy_gem_good_vertex[gemIndex]->Fill(xz_separation[gemIndex], y_separation);
          h_dx_gemX[gemIndex]->Fill(xz_separation[gemIndex], trk_x_local[gemIndex][j]);
          h_gem_xz_good_vertex[gemIndex]->Fill(p_gem_hit.X(), p_gem_hit.Z());
          h_intersection_xz_good_vertex[gemIndex]->Fill(interaction.X(), interaction.Z());
        } else {
          h_dx_dy_gem_bad_vertex[gemIndex]->Fill(xz_separation[gemIndex], y_separation);
          // h_dx_gemX[gemIndex]->Fill(xz_separation[gemIndex], p_gem_hit.X());
          h_gem_xz_bad_vertex[gemIndex]->Fill(p_gem_hit.X(), p_gem_hit.Z());
          h_intersection_xz_bad_vertex[gemIndex]->Fill(interaction.X(), interaction.Z());
        }
        h_intersection_xz[gemIndex]->Fill(interaction.X(), interaction.Z());
        h_gem1dx_vs_targetz->Fill(xz_separation[0], p_react.Z());
      }

      h_gem1dx_vs_gem2dx->Fill(xz_separation[0], xz_separation[1]);
      h_gem1dx_vs_plane->Fill(hodoPlane[j], xz_separation[0]);
      h_gem1dx_vs_paddle[hodoPlane[j]]->Fill(hodoPaddle[j], xz_separation[0]);
      h_hodo_hit_xz->Fill(p_hit.X(), p_hit.Z());
      h_target_xz->Fill(p_react.X(), p_react.Z());
      if (good_vertez) {
        h_gem1dx_vs_gem2dx_good->Fill(xz_separation[0], xz_separation[1]);
      } else {
        h_gem1dx_vs_gem2dx_bad->Fill(xz_separation[0], xz_separation[1]);
      }
    }
    // Print the status as a percentage
    if (i % (nEntries / 100) == 0) {
      std::cout << "\rProcessing: " << int(i * 100.0 / nEntries) << "% completed." << std::flush;
    }
  }
  // End Event Loop
  std::cout << "\nProcessing completed." << std::endl;
  // Close the file

  // Create output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return;
  }

  // Create projections of dx_dy (all, good, and bad) onto the x and y axes
  TH1D *h_dx_dy_proj_x[2], *h_dx_dy_proj_y[2];
  TH1D *h_dx_dy_good_proj_x[2], *h_dx_dy_good_proj_y[2];
  TH1D *h_dx_dy_bad_proj_x[2], *h_dx_dy_bad_proj_y[2];

  for (int gemIndex = 0; gemIndex < 2; ++gemIndex) {
    h_dx_dy_proj_x[gemIndex] = h_dx_dy_gem[gemIndex]->ProjectionX(Form("h_dx_dy_proj_x_gem%d", gemIndex));
    h_dx_dy_proj_x[gemIndex]->SetTitle(Form("Projection of dx_dy GEM %d on X;dx (cm);Counts", gemIndex));
    h_dx_dy_proj_y[gemIndex] = h_dx_dy_gem[gemIndex]->ProjectionY(Form("h_dx_dy_proj_y_gem%d", gemIndex));
    h_dx_dy_proj_y[gemIndex]->SetTitle(Form("Projection of dx_dy GEM %d on Y;dy (cm);Counts", gemIndex));

    h_dx_dy_good_proj_x[gemIndex] =
        h_dx_dy_gem_good_vertex[gemIndex]->ProjectionX(Form("h_dx_dy_good_proj_x_gem%d", gemIndex));
    h_dx_dy_good_proj_x[gemIndex]->SetTitle(
        Form("Projection of dx_dy (Good Vertex) GEM %d on X;dx (cm);Counts", gemIndex));
    h_dx_dy_good_proj_y[gemIndex] =
        h_dx_dy_gem_good_vertex[gemIndex]->ProjectionY(Form("h_dx_dy_good_proj_y_gem%d", gemIndex));
    h_dx_dy_good_proj_y[gemIndex]->SetTitle(
        Form("Projection of dx_dy (Good Vertex) GEM %d on Y;dy (cm);Counts", gemIndex));

    h_dx_dy_bad_proj_x[gemIndex] =
        h_dx_dy_gem_bad_vertex[gemIndex]->ProjectionX(Form("h_dx_dy_bad_proj_x_gem%d", gemIndex));
    h_dx_dy_bad_proj_x[gemIndex]->SetTitle(
        Form("Projection of dx_dy (Bad Vertex) GEM %d on X;dx (cm);Counts", gemIndex));
    h_dx_dy_bad_proj_y[gemIndex] =
        h_dx_dy_gem_bad_vertex[gemIndex]->ProjectionY(Form("h_dx_dy_bad_proj_y_gem%d", gemIndex));
    h_dx_dy_bad_proj_y[gemIndex]->SetTitle(
        Form("Projection of dx_dy (Bad Vertex) GEM %d on Y;dy (cm);Counts", gemIndex));
  }
  // Write the projections to the output file
  outputFile->mkdir("dx_dy");
  outputFile->cd("dx_dy");
  for (int gemIndex = 0; gemIndex < 2; ++gemIndex) {
    h_dx_dy_proj_x[gemIndex]->Write();
    h_dx_dy_proj_y[gemIndex]->Write();
    h_dx_dy_good_proj_x[gemIndex]->Write();
    h_dx_dy_good_proj_y[gemIndex]->Write();
    h_dx_dy_bad_proj_x[gemIndex]->Write();
    h_dx_dy_bad_proj_y[gemIndex]->Write();
    h_dx_dy_gem[gemIndex]->Write();
    h_dx_dy_gem_good_vertex[gemIndex]->Write();
    h_dx_dy_gem_bad_vertex[gemIndex]->Write();
    h_gem_xz[gemIndex]->Write();
    h_gem_xz_good_vertex[gemIndex]->Write();
    h_gem_xz_bad_vertex[gemIndex]->Write();
    h_intersection_xz[gemIndex]->Write();
    h_intersection_xz_good_vertex[gemIndex]->Write();
    h_intersection_xz_bad_vertex[gemIndex]->Write();
    h_dx_gemX[gemIndex]->Write();
  }
  h_gem1dx_vs_plane->Write();
  for (int planeIndex = 0; planeIndex < nPlanes; ++planeIndex) {
    h_gem1dx_vs_paddle[planeIndex]->Write();
  }

  h_gem1dx_vs_gem2dx->Write();
  h_gem1dx_vs_gem2dx_good->Write();
  h_gem1dx_vs_gem2dx_bad->Write();
  h_gem1dx_vs_targetz->Write();
  h_trk_projz->Write();
  h_projz_vs_vertexz->Write();
  // Write the hodoscope hit xz histogram
  outputFile->mkdir("diagnostic_plots");
  outputFile->cd("diagnostic_plots");
  h_hodo_hit_xz->Write();
  h_target_xz->Write();
  // // Close the output file
  outputFile->Close();
  file->Close();
  // delete outputFile;
  // // Clean up
  // delete h_dx_dy;
  // delete h_dx_dy_good_vertex;
  // delete h_dx_dy_bad_vertex;
}