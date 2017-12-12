<?php
  $suffix = $_GET["pic"];
  $filename = "/tmp/MUT" . $suffix;
  header('Content-Type: image/png');
  $img = imagecreatefrompng($filename);
  imagepng($img);
  imagedestroy($img);
?>
