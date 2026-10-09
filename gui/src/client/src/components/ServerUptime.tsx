import React from "react";
import Stack from "@mui/material/Stack";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import InfoOutlinedIcon from "@mui/icons-material/InfoOutlined";
import { useQuery } from "@tanstack/react-query";
import { getServerStartup } from "../features/server/api";
import { formatUptime } from "../features/server/utils";

// The startup epoch this reads (see /api/startup) is stored in Redis, the
// same store sessions live in -- it only resets when Redis itself was
// restarted/flushed, not on every backend process restart. So a short uptime
// here is a real signal that older sessions may be gone, not just noise.
export function ServerUptime() {
  const { data, isLoading, isError } = useQuery({
    queryKey: ["serverStartup"],
    queryFn: ({ signal }) => getServerStartup(signal),
    refetchInterval: 60_000,
    retry: false,
  });

  if (isLoading || isError || !data) return null;

  return (
    <Tooltip
      title="If this looks shorter than expected, the server was recently restarted and older sessions may no longer be retrievable."
      arrow
    >
      <Stack
        direction="row"
        spacing={0.5}
        alignItems="center"
        sx={{ cursor: "help", color: "text.secondary", width: "fit-content" }}
      >
        <InfoOutlinedIcon fontSize="inherit" />
        <Typography variant="caption">
          Server up for {formatUptime(data.uptime)}
        </Typography>
      </Stack>
    </Tooltip>
  );
}
