import React from "react";
import type { NotificationSeverity } from "../features/notifications/types";

// Extend the Notification interface with a "level" property
export interface Notification {
  id: number;
  content: string;
  timestamp: Date;
  isNew: boolean;
  level: NotificationSeverity;
}

export interface NotificationContextType {
  notifications: Notification[];
  pushNotification: (content: string, level?: NotificationSeverity) => void;
  markAllAsRead: () => void;
  clearNotifications: () => void;
}

const NotificationContext = React.createContext<NotificationContextType>({
  notifications: [],
  pushNotification: () => {},
  markAllAsRead: () => {},
  clearNotifications: () => {},
})

export function useNotifications() {
  return React.useContext(NotificationContext);
}

export const NotificationProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [notifications, setNotifications] = React.useState<Notification[]>([]);

  const pushNotification = React.useCallback(
    (content: string, level: NotificationSeverity = "info") => {
      console.log("Pushing notification:", content, level);
      setNotifications((prev) => {
        const newNotification: Notification = {
          id: Date.now(),
          content,
          timestamp: new Date(),
          isNew: true,
          level,
        };

        let updated = [...prev, newNotification];

        if (updated.length > 50) {
          updated = updated.slice(updated.length - 50);
        }
        return updated;
      })
    },
    []
  )

  const markAllAsRead = React.useCallback(() => {
    setNotifications((prev) => {
      const hasUnread = prev.some((notification) => notification.isNew);
      if (!hasUnread) return prev;
      return prev.map((notification) => ({ ...notification, isNew: false }))
    });
  }, [])

  const clearNotifications = React.useCallback(() => {
    setNotifications([]);
  }, [])

  const contextValue = React.useMemo(
    () => ({ notifications, pushNotification, markAllAsRead, clearNotifications }),
    [notifications, pushNotification, markAllAsRead, clearNotifications]
  )

  return (
    <NotificationContext.Provider value={contextValue}>
      {children}
    </NotificationContext.Provider>
  )
}