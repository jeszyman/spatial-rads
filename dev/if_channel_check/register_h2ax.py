"""
Register H2AX IHC (left panel of slide 5) to CosMx IF composite (right panel)
using tumor outline matching. Output a rectified H2AX PNG in CosMx coord frame.
"""
import numpy as np
from PIL import Image
from skimage import color, filters, measure, morphology, transform, io
import matplotlib.pyplot as plt

SLIDE = "/tmp/h2ax_register/slide5-5.png"
OUT = "/tmp/h2ax_register"

# ---- Step 1: crop left (H2AX) and right (CosMx IF) panels ----
slide = np.array(Image.open(SLIDE).convert("RGB"))
H, W = slide.shape[:2]
# Heuristic split: left panel in left half, right panel in right half, below title
y0 = int(0.18 * H)
y1 = int(0.95 * H)
left = slide[y0:y1, 0:int(0.50 * W)]
right = slide[y0:y1, int(0.50 * W):W]
Image.fromarray(left).save(f"{OUT}/panel_h2ax_raw.png")
Image.fromarray(right).save(f"{OUT}/panel_cosmx_raw.png")
print(f"Left panel: {left.shape}, right panel: {right.shape}")

# ---- Step 2: segment tumor silhouette in each panel ----
def segment_tumor_h2ax(img_rgb):
    """H2AX panel: tumor is pink/brown on white background; segment by color distance from white."""
    # Distance from white in RGB
    dist_from_white = np.sqrt(np.sum((img_rgb.astype(float) - 255.0) ** 2, axis=-1))
    # Tissue has dist > some threshold (pinkish-brown differs from pure white)
    thr = 30  # gentle threshold
    mask = dist_from_white > thr
    mask = morphology.remove_small_holes(mask, area_threshold=20000)
    mask = morphology.remove_small_objects(mask, min_size=50000)
    lbl = measure.label(mask)
    if lbl.max() == 0:
        return mask
    regs = measure.regionprops(lbl)
    biggest = max(regs, key=lambda r: r.area)
    return (lbl == biggest.label)

def segment_tumor_cosmx(img_rgb):
    """CosMx panel: tumor is bright on black background."""
    gray = color.rgb2gray(img_rgb)
    thr = max(filters.threshold_otsu(gray), 0.05)
    mask = gray > thr
    mask = morphology.remove_small_holes(mask, area_threshold=10000)
    mask = morphology.remove_small_objects(mask, min_size=50000)
    lbl = measure.label(mask)
    if lbl.max() == 0:
        return mask
    regs = measure.regionprops(lbl)
    biggest = max(regs, key=lambda r: r.area)
    return (lbl == biggest.label)

mask_h = segment_tumor_h2ax(left)
mask_c = segment_tumor_cosmx(right)

Image.fromarray((mask_h * 255).astype(np.uint8)).save(f"{OUT}/mask_h2ax.png")
Image.fromarray((mask_c * 255).astype(np.uint8)).save(f"{OUT}/mask_cosmx.png")
print(f"H2AX mask area: {mask_h.sum()} px, CosMx mask area: {mask_c.sum()} px")

# ---- Step 3: extract ordered contour points from each mask ----
def ordered_contour(mask, n_points=200):
    contours = measure.find_contours(mask.astype(float), level=0.5)
    contour = max(contours, key=len)  # biggest contour
    # Subsample evenly
    idx = np.linspace(0, len(contour) - 1, n_points).astype(int)
    return contour[idx]  # shape (n_points, 2), (row, col) i.e. (y, x)

pts_h = ordered_contour(mask_h, n_points=300)
pts_c = ordered_contour(mask_c, n_points=300)
print(f"H2AX contour pts: {pts_h.shape}, CosMx contour pts: {pts_c.shape}")

# ---- Step 4: align via shape-centered rotation search + affine fit ----
# Strategy: center both contours on their centroid, then rotate H2AX through
# 360 degrees in 1-deg steps, for each rotation compute the sum-of-min-distance
# between H2AX points and CosMx points. Pick best rotation.
def center_and_scale(pts):
    c = pts.mean(axis=0)
    d = pts - c
    scale = np.sqrt((d ** 2).sum(axis=1)).mean()
    return d / scale, c, scale

ph, ch, sh = center_and_scale(pts_h)
pc, cc, sc = center_and_scale(pts_c)

def rotate_pts(pts, theta):
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    # pts are (y, x). Rotate in xy -> swap
    xy = pts[:, [1, 0]]
    xy_r = xy @ R.T
    return xy_r[:, [1, 0]]

def nn_loss(a, b):
    # For each point in a, distance to nearest in b
    d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)
    return np.sqrt(d2.min(1)).mean() + np.sqrt(d2.min(0)).mean()

best = (np.inf, 0)
for deg in range(-180, 180, 1):
    theta = np.deg2rad(deg)
    ph_r = rotate_pts(ph, theta)
    loss = nn_loss(ph_r, pc)
    if loss < best[0]:
        best = (loss, deg)
print(f"Best rotation: {best[1]} deg, loss {best[0]:.4f}")

# Refine at 0.1 deg resolution
best2 = (best[0], best[1])
for deg10 in range(-20, 21):
    deg = best[1] + deg10 * 0.1
    theta = np.deg2rad(deg)
    ph_r = rotate_pts(ph, theta)
    loss = nn_loss(ph_r, pc)
    if loss < best2[0]:
        best2 = (loss, deg)
theta_fit = np.deg2rad(best2[1])
print(f"Refined rotation: {best2[1]} deg, loss {best2[0]:.4f}")

# ---- Step 5: apply rigid transform to raw H2AX panel ----
# Build the transform: translate H2AX centroid to origin, rotate by theta, scale by sc/sh, translate to CosMx centroid
s_ratio = sc / sh
# skimage.transform.AffineTransform operates in (x, y) order
# forward mapping: xy_cosmx = (xy_h2ax - ch_xy) * R * s + cc_xy
ch_xy = ch[[1, 0]]
cc_xy = cc[[1, 0]]

# skimage transforms operate on images; we want to warp the H2AX raw panel into CosMx frame
# Using inverse_map: for each target pixel (x, y) in cosmx frame, compute source pixel in h2ax frame
from skimage.transform import AffineTransform, warp
# forward: cosmx = R*s*(h2ax - ch) + cc
# so inverse: h2ax = R^-1 * (cosmx - cc) / s + ch
inv_R = np.array([[np.cos(-theta_fit), -np.sin(-theta_fit)],
                  [np.sin(-theta_fit),  np.cos(-theta_fit)]])
def inv_map(xy):
    xy_rel = xy - cc_xy[None, :]
    xy_rot = xy_rel @ inv_R.T / s_ratio
    return xy_rot + ch_xy[None, :]

target_shape = right.shape[:2]
warped = warp(left, inverse_map=inv_map, output_shape=target_shape, order=1, mode="constant", cval=1.0)
warped_u8 = (warped * 255).astype(np.uint8) if warped.dtype == float else warped.astype(np.uint8)
Image.fromarray(warped_u8).save(f"{OUT}/h2ax_rectified.png")
print(f"Saved rectified H2AX: {OUT}/h2ax_rectified.png shape {warped_u8.shape}")

# ---- Step 6: diagnostic overlay: original CosMx IF + rectified H2AX side by side and overlaid ----
fig, ax = plt.subplots(1, 3, figsize=(18, 6))
ax[0].imshow(left); ax[0].set_title("H2AX raw"); ax[0].axis("off")
ax[1].imshow(right); ax[1].set_title("CosMx IF composite"); ax[1].axis("off")
ax[2].imshow(right); ax[2].imshow(warped_u8, alpha=0.5); ax[2].set_title(f"Rectified H2AX over CosMx ({best2[1]:.1f} deg)"); ax[2].axis("off")
plt.tight_layout()
plt.savefig(f"{OUT}/registration_diagnostic.png", dpi=150)
print(f"Saved diagnostic: {OUT}/registration_diagnostic.png")

# Save transform params for reuse
np.savez(f"{OUT}/h2ax_transform.npz",
         theta_rad=theta_fit, scale=s_ratio,
         centroid_h2ax=ch_xy, centroid_cosmx=cc_xy,
         h2ax_panel_shape=left.shape, cosmx_panel_shape=right.shape)
print(f"Saved transform: {OUT}/h2ax_transform.npz")
