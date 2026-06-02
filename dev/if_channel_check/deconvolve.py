from skimage import io, color, exposure
import numpy as np

img = io.imread("/tmp/mbrt-4hr.png")
if img.shape[-1] == 4:
    img = img[..., :3]
hed = color.rgb2hed(img)
dab = hed[..., 2]
dab_norm = exposure.rescale_intensity(dab, out_range=(0, 1))
p1, p99 = np.percentile(dab_norm, (1, 99))
dab_stretched = exposure.rescale_intensity(dab_norm, in_range=(p1, p99))

# Save as grayscale where darker = higher H2AX (matches visual intuition of stain intensity)
out_gray = (dab_stretched * 255).astype(np.uint8)
io.imsave("/tmp/mbrt-4hr-dab.png", out_gray)

# Also save an inverted high-contrast version where bright = high H2AX
out_hot = (dab_stretched * 255).astype(np.uint8)
io.imsave("/tmp/mbrt-4hr-dab-bright.png", out_hot)

print(f"DAB image: {dab_stretched.shape}, range [{dab_stretched.min():.3f}, {dab_stretched.max():.3f}]")
print("Saved /tmp/mbrt-4hr-dab.png and /tmp/mbrt-4hr-dab-bright.png")
