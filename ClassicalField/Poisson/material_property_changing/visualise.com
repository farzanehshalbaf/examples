#Read in the sequence of nodal positions.
for $i (0..2)
  {
	 $filename = sprintf("./output/TIMESTP_%01d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elem ./output/TIMESTP_0.part0.exelem;

gfx define faces egroup PoissonRegion
gfx create window 1
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range -5.0 5.0 extend_above extend_below rainbow colour_range 0 1 component 1
gfx modify spectrum default linear reverse range -5. 5.0 extend_above extend_below banded number_of_bands 20 band_ratio 0.06 component 1
gfx modify g_element PoissonRegion node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 
gfx modify g_element PoissonRegion lines select_on material default selected_material default_selected
gfx modify g_element PoissonRegion surfaces select_on material default data Phi spectrum default selected_material default_selected render_shaded




gfx timekeeper default set 1.0;
gfx timekeeper default play;
gfx create time_editor
