import React from "react";
import Chip from "@mui/material/Chip";

interface ItemKindChipProps {
  itemKind: string;
}

export const ItemKindChip: React.FC<ItemKindChipProps> = ({ itemKind }) => {
  const labelMap: Record<string, string> = {
    "compound": "Compound",
    "gene_cluster": "BGC",
  };

  const colorMap: Record<string, "default" | "primary" | "secondary"> = {
    "compound": "primary",
    "gene_cluster": "secondary",
  };

  return (
    <Chip
      label={labelMap[itemKind] || "Unknown"}
      color={colorMap[itemKind] || "default"}
      size="small"
      sx={{ fontSize: "0.7rem", height: 20 }}
    />
  );
}