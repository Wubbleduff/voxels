





// Goal:
// Generate terrain with view distance 262,144 m
//
// Need a generated region of 524288^3 = 1.44115188075855872 * 10^17
// To simplify terrain, will reduce LOD as distance increase from player.
//
// Example:
// LOD0_WIDTH = 8
// MAX_LOD = 3
//
//      *               *               *               *               *               *               *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *               *               *               *               *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *       *       *       *       *       *       *       *       *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                      *       *       *       *       *       *       *       *       *                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *       *       *   *   *   *   *   *   *   *   *       *       *               *               *     
//                                                                      |                                                                     
//                                                      *   *   *   *   *   *   *   *   *                                                     
//                                                                      |                                                                     
//                                      *       *       *   *   * * * * * * * * *   *   *       *       *                                     
//                                                              * * * * * * * * *                                                             
//                                                      *   *   * * * * * * * * *   *   *                                                     
//                                                              * * * * * * * * *                                                             
//  ----*---------------*---------------*-------*-------*---*---*-*-*-*-*-*-*-*-*---*---*-------*-------*---------------*---------------*-------------------
//                                                              * * * * * * * * *                                                             
//                                                      *   *   * * * * * * * * *   *   *                                                     
//                                                              * * * * * * * * *                                                             
//                                      *       *       *   *   * * * * * * * * *   *   *       *       *                                     
//                                                                      |                                                                     
//                                                      *   *   *   *   *   *   *   *   *                                                     
//                                                                      |                                                                     
//      *               *               *       *       *   *   *   *   *   *   *   *   *       *       *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                      *       *       *       *       *       *       *       *       *                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *       *       *       *       *       *       *       *       *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *               *               *               *               *               *               *     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//                                                                      |                                                                     
//      *               *               *               *               *               *               *               *               *     
//                                                    9 8 7 6 5 4 3 2 1 | 1 2 3 4 5 6 7 8 9                                                   
// 
// In 3D:
// Each LOD is 2x wider / sparser than the last
// (just store all samples in each layer for simplicity)
// # samples = 8^3 + 8^3 + 8^3 + 8^3 = 2048
//
// Example:
// LOD0_width = 64
// MAX_LOD = 18
// # samples = 18 * 64^3 = 4718592 (~5 million)
//

static const u32 LOD0_WIDTH = 64;
static const u32 MAX_LOD = 18;

struct TerrainSamples
{
    float samples[LOD0_WIDTH*LOD0_WIDTH*LOD0_WIDTH];
}

// Input:
// * LOD 0 width (static u32)
// * Max LOD (static u32)
// * Camera pos (3x s32)
// Output:
// * List of tris (vertices, normals, indices)
// * TODO Acceleration spatial lookup structure for finding terrain
void generate_terrain(
        f32_m* out_vx,
        f32_m* out_vy,
        f32_m* out_vz,
        f32_m* out_nx,
        f32_m* out_ny,
        f32_m* out_nz,
        u32_m* out_indices,
        s32 cam_x,
        s32 cam_y,
        s32 cam_z)
{
}




