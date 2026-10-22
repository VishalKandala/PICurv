- `picurv validate` warns when a post recipe sets `io.output_directory`, which is ignored: output is filed by
  recipe identity, so recipes differing only in it overwrote each other. Shipped recipes no longer set it.
